function [X,Y,C] = updateXYC2(X, Y, C, GwX, GwY, GwC, GqX, GqY, GqC, step_R, step_p, nData)

minstep_rot = 1e-5;
minstep_trans = 1e-5;

dwX = step_R/nData*GwX;
dwX = dwX * (1 + minstep_rot/(norm(dwX)+1e-9));

dqX = step_p/nData*GqX';
dqX = dqX * (1 + minstep_trans/(norm(dqX)+1e-9));

dwY = step_R/nData*GwY;
dwY = dwY * (1 + minstep_rot/(norm(dwY)+1e-9));

dqY = step_p/nData*GqY';
dqY = dqY * (1 + minstep_trans/(norm(dqY)+1e-9));

X = X * [LargeSO3(dwX), dqX; 0,0,0,1];
Y = Y * [LargeSO3(dwY), dqY; 0,0,0,1];

n = size(C,3);
for i = 1:n
    dwC = step_R*GwC(i,:);
    dwC = dwC * (1 + minstep_rot/(norm(dwC)+1e-9));
    
    dqC = step_p*GqC(i,:)';
    dqC = dqC * (1 + minstep_trans/(norm(dqC)+1e-9));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     C{i} = C{i} * [LargeSO3(dwC), dqC; 0,0,0,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    w = dwC;
    TINY=1e-14;
    theta = norm(w);
	if ( abs(theta) < TINY ) 
        continue; 
    end
    w=w/theta;
	st = sin(theta);
    ct = cos(theta);
    vt = 1.0 - ct;
    t0 = w(3) * st;
    t1 = w(2) * st;
    t2 = w(1) * st;
    w0=w(1);w1=w(2);w2=w(3);

    T = [w0 * w0 * vt + ct, w0 * w1 * vt - t0, w0 * w2 * vt + t1, dqC(1);
        w0 * w1 * vt + t0, w1 * w1 * vt + ct, w1 * w2 * vt - t2, dqC(2);
        w0 * w2 * vt - t1, w1 * w2 * vt + t2, w2 * w2 * vt + ct, dqC(3);
        0, 0, 0, 1];
    
    C{i} = C{i} * T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end


end

