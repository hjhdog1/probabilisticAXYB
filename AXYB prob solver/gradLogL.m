function [GwX, GwY, GwC, GqX, GqY, GqC, J] = gradLogL(X, Y, A, B, C, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf)
% noiseConf == 1 for N = CX^{-1}A^{-1} and M = C^{-1}YB
% noiseConf == 2 for N = XC^{-1}A and M = C^{-1}YB

n = size(A,3);
if n == 1
    A = {A};
    B = {B};
    C = {C};
    invSig_wN = {invSig_wN};
    invSig_pN = {invSig_pN};
    invSig_wM = {invSig_wM};
    invSig_pM = {invSig_pM};
end

RX = X(1:3,1:3);
RY = Y(1:3,1:3);
pX = X(1:3,4);
pY = Y(1:3,4);
invX = invSE3(X);

N = cell(1,n);
M = cell(1,n);

GwX = zeros(1,3);
GwY = zeros(1,3);
GwC = zeros(n,3);
GqX = zeros(1,3);
GqY = zeros(1,3);
GqC = zeros(n,3);
J = 0;
for i = 1:n

    RA = A{i}(1:3,1:3);
%     RB = B{i}(1:3,1:3);
    RC = C{i}(1:3,1:3);
    pA = A{i}(1:3,4);
    pB = B{i}(1:3,4);
    pC = C{i}(1:3,4);
    
    RC_trans = RC';
    invC = [RC_trans, -RC_trans*pC; 0,0,0,1];
    if noiseConf == 1
        N{i} = C{i}*invX*invSE3(A{i});
%         N{i} = A{i}*X*invC;
    elseif noiseConf == 2
        N{i} = X*invC*A{i};
    else
        error('Error: check the noise configuration.')
    end        
%     wN = LogSO3(N{i}(1:3,1:3));
    wN = LogSO3(N{i});
    pN = N{i}(1:3,4);

    M{i} = invC*Y*B{i};
%     wM = LogSO3(M{i}(1:3,1:3));
    wM = LogSO3(M{i});
    pM = M{i}(1:3,4);
    
    WN = invSig_wN{i};
    PN = invSig_pN{i};
    WM = invSig_wM{i};
    PM = invSig_pM{i};
    %%%%%%%%%%%%%%%%%%%%%%%

    if noiseConf == 1
        invD1 = diffExpInv_rotation_only(wN);
        invD2 = diffExpInv_rotation_only(wM);
        
        u = pC-pN;
        U = -wN'*WN*invD1*RC - pN'*PN*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0]*RC;
        u = RC'*(RY*pB);
        V = (wM'*WM)*invD2 - pM'*PM*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
        
        GwX = GwX - U;
        GwY = GwY - V*RC'*RY;        
        u = RC'*(pY-pC);
        GwC(i,:) = U + V - pM'*PM*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
        GqX = GqX + pN'*PN*RC;
        GqY = GqY - pM'*PM*RC'*RY;
        GqC(i,:) = -pN'*PN*RC + pM'*PM;
        
    elseif noiseConf == 2
        invD1 = diffExpInv_rotation_only(wN);
        invD2 = diffExpInv_rotation_only(wM);
        
        u = RX*RC'*(pA-pC);
        S = (wN'*WN)*invD1*RX - pN'*PN*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0]*RX;
        u = RC'*RY*pB;
        V = (wM'*WM)*invD2 - pM'*PM*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
        
        GwX = GwX - S;
        GwY = GwY - V*RC'*RY;        
        u = RC'*(pY-pC);
        GwC(i,:) = S + V - pM'*PM*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
        GqX = GqX - pN'*PN*RX;
        GqY = GqY - pM'*PM*RC'*RY;
        GqC(i,:) = pN'*PN*RX + pM'*PM;
    end
    J = J - 0.5*(wN'*WN*wN + pN'*PN*pN + wM'*WM*wM + pM'*PM*pM);
    
end

end

