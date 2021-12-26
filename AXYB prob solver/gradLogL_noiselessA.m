function [GwX, GwY, GqX, GqY, J] = gradLogL_noiselessA(X, Y, A, B, invSig_wM, invSig_pM)
% M = X^{-1}A^{-1}YB

n = size(A,3);
if n == 1
    A = {A};
    B = {B};
    invSig_wM = {invSig_wM};
    invSig_pM = {invSig_pM};
end

RX = X(1:3,1:3);
RY = Y(1:3,1:3);
pX = X(1:3,4);
pY = Y(1:3,4);
invX = invSE3(X);

M = cell(1,n);

GwX = zeros(1,3);
GwY = zeros(1,3);
GqX = zeros(1,3);
GqY = zeros(1,3);
J = 0;
for i = 1:n

    RA = A{i}(1:3,1:3);
    RB = B{i}(1:3,1:3);
    pA = A{i}(1:3,4);
    pB = B{i}(1:3,4);
    
    M{i} = invX*invSE3(A{i})*Y*B{i};
    RM = M{i}(1:3,1:3);
    wM = LogSO3(RM);
    pM = M{i}(1:3,4);
    
    WM = invSig_wM{i};
    PM = invSig_pM{i};
    
    % computation of gradients
    invD = diffExpInv_rotation_only(wM);
    H = wM'*WM*invD;       % D_i in paper
    pMPM = pM'*PM;
    
    u = pM;
    GwX = GwX + H - pMPM*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
    
    u = pB;
    GwY = GwY - H*RM*RB' + pMPM*RM*RB'*[0, -u(3), u(2); u(3), 0, -u(1); -u(2), u(1), 0];
    
    GqX = GqX + pMPM;
    GqY = GqY - pMPM*RM*RB';
        
    % computation of cost function
    J = J - 0.5*(wM'*WM*wM + pM'*PM*pM);
    
end

end

