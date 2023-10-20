function R = LargeSO3(w)

TINY=1e-14;
if(size(w)==[3,1] | size(w)==[1,3])
    theta = norm(w);
	if ( abs(theta) < TINY ) 
        R=eye(3);
        return 
    end
    w=w/theta;
	st = sin(theta);
    ct = cos(theta);
    vt = 1.0 - ct;
    t0 = w(3) * st;
    t1 = w(2) * st;
    t2 = w(1) * st;
    w0=w(1);w1=w(2);w2=w(3);

    R = [w0 * w0 * vt + ct, w0 * w1 * vt - t0, w0 * w2 * vt + t1;
        w0 * w1 * vt + t0, w1 * w1 * vt + ct, w1 * w2 * vt - t2;
        w0 * w2 * vt - t1, w1 * w2 * vt + t2, w2 * w2 * vt + ct];

% 
%     R = zeros(3);
%     R(1,1)=w0 * w0 * vt + ct;
%     R(2,1)=w0 * w1 * vt + t0;
%     R(3,1)=w0 * w2 * vt - t1;
%     R(1,2)=w0 * w1 * vt - t0;
%     R(2,2)=w1 * w1 * vt + ct;
%     R(3,2)=w1 * w2 * vt + t2;
%     R(1,3)=w0 * w2 * vt + t1;
%     R(2,3)=w1 * w2 * vt - t2;
%     R(3,3)=w2 * w2 * vt + ct;
else
    error('SO3(w) w is not 3 x 1')
end