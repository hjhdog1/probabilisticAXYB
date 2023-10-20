function [covX,covY, Z] = computeUncertainty_noiseConf3(X, Y, B, invSig_wM, invSig_pM)



n = size(B,3);

%% Q
Qmx = zeros(6*n,6);
Qmy = zeros(6*n,6);

for i = 1:n
    RB = B(1:3,1:3,i);
    pB = B(1:3,4,i);
    
    idx = 6*(i-1)+1:6*(i-1)+3;
    
    
    Qmx(idx, 1:3) = -eye(3);
    Qmx(idx+3, 4:6) = -eye(3);
    
    
    Qmy(idx, 1:3) = RB';
    Qmy(idx+3, 1:3) = -RB' * skew(pB);
    Qmy(idx+3, 4:6) = RB';
    
end

Q = [Qmx, Qmy];

%% Z
E_true = eye(6*n);
for i = 1:n    
    
    idx = 6*(i-1)+1:6*(i-1)+3;
    
    E_true(idx,idx) = inv(invSig_wM(:,:,i));
    E_true(idx+3,idx+3) = inv(invSig_pM(:,:,i));
end

invE = eye(6*n);
for i = 1:n    
    
    idx = 6*(i-1)+1:6*(i-1)+3;
    
    invE(idx,idx) = invSig_wM(:,:,i);
    invE(idx+3,idx+3) = invSig_pM(:,:,i);
    
end

Z = -(Q'*invE*Q)\Q'*invE;

E_XYC = Z * E_true * Z';
covX = E_XYC(1:6,1:6);
covY = E_XYC(7:12,7:12);


end

