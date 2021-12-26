function J = computeLogL_noiselessA(X, Y, A, B, invSig_wM, invSig_pM)
% M = X^{-1}A^{-1}YB

invX = invSE3(X);

n = size(A,3);
M = zeros(4,4,n);

J = 0;
for i = 1:n
    
    
    M(:,:,i) = invX*invSE3(A(:,:,i))*Y*B(:,:,i);
    wM = LogSO3(M(1:3,1:3,i));
    pM = M(1:3,4,i);
    
    WM = invSig_wM(:,:,i);
    PM = invSig_pM(:,:,i);
        
    J = J - 0.5*(wM'*WM*wM + pM'*PM*pM);
    
end

end

