function [covX,covY, Z] = computeUncertainty(X, Y, C, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf)



n = size(C,3);

%% Q
Qnx = zeros(6*n,6);
Qny = zeros(6*n,6);
Qnc = zeros(6*n,6*n);

Qmx = zeros(6*n,6);
Qmy = zeros(6*n,6);
Qmc = zeros(6*n,6*n);

RX = X(1:3,1:3);
pX = X(1:3,4);
RY = Y(1:3,1:3);
pY = Y(1:3,4);

for i = 1:n
    RC = C(1:3,1:3,i);
    pC = C(1:3,4,i);
    
    
    idx = 6*(i-1)+1:6*(i-1)+3;
    
    if noiseConf == 1
        Qnx(idx, 1:3) = -RC;
        Qnx(idx+3, 1:3) = -skew(pC)*RC;
        Qnx(idx+3, 4:6) = -RC;

        Qnc(idx, idx) = RC;
        Qnc(idx+3, idx) = skew(pC)*RC;
        Qnc(idx+3, idx+3) = RC;
    elseif noiseConf == 2
        Qnx(idx, 1:3) = RX;
        Qnx(idx+3, 1:3) = skew(pX)*RX;
        Qnx(idx+3, 4:6) = RX;

        Qnc(idx, idx) = -RX;
        Qnc(idx+3, idx) = -skew(pX)*RX;
        Qnc(idx+3, idx+3) = -RX;
    else
        error('Check noise configuration.')
    end
    
    Qmy(idx, 1:3) = RC'*RY;
    Qmy(idx+3, 1:3) = skew(RC'*(pY-pC))*RC'*RY;
    Qmy(idx+3, 4:6) = RC'*RY;
    
    Qmc(idx, idx) = -eye(3);
    Qmc(idx+3, idx+3) = -eye(3);
    
end

Q = [Qnx, Qny, Qnc;
     Qmx, Qmy, Qmc];

%% Z
E_true = eye(12*n);
for i = 1:n    
    
    idx = 6*(i-1)+1:6*(i-1)+3;
    
    E_true(idx,idx) = inv(invSig_wN(:,:,i));
    E_true(idx+3,idx+3) = inv(invSig_pN(:,:,i));
    E_true(idx + 6*n, idx + 6*n) = inv(invSig_wM(:,:,i));
    E_true(idx+3 + 6*n, idx+3 + 6*n) = inv(invSig_pM(:,:,i));
end

invE = eye(12*n);
for i = 1:n    
    
    idx = 6*(i-1)+1:6*(i-1)+3;
    
    invE(idx,idx) = invSig_wN(:,:,i);
    invE(idx+3,idx+3) = invSig_pN(:,:,i);
    invE(idx + 6*n, idx + 6*n) = invSig_wM(:,:,i);
    invE(idx+3 + 6*n, idx+3 + 6*n) = invSig_pM(:,:,i);
    
end

Z = -(Q'*invE*Q)\Q'*invE;

E_XYC = Z * E_true * Z';
covX = E_XYC(1:6,1:6);
covY = E_XYC(7:12,7:12);


end

