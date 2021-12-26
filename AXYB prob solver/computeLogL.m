function J = computeLogL(X, Y, A, B, C, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf)
% noiseConf == 1 for N = CX^{-1}A^{-1} and M = C^{-1}YB
% noiseConf == 2 for N = XC^{-1}A and M = C^{-1}YB

n = size(A,3);

N = zeros(4,4,n);
M = zeros(4,4,n);

J = 0;
for i = 1:n
    
    if noiseConf == 1
        N(:,:,i) = C(:,:,i)*invSE3(X)*invSE3(A(:,:,i));
    elseif noiseConf == 2
        N(:,:,i) = X*invSE3(C(:,:,i))*A(:,:,i);
    else
        error('Error: check the noise configuration.')
    end        
    wN = LogSO3(N(1:3,1:3,i));
    pN = N(1:3,4,i);

    M(:,:,i) = invSE3(C(:,:,i))*Y*B(:,:,i);
    wM = LogSO3(M(1:3,1:3,i));
    pM = M(1:3,4,i);
    
    WN = invSig_wN(:,:,i);
    PN = invSig_pN(:,:,i);
    WM = invSig_wM(:,:,i);
    PM = invSig_pM(:,:,i);
        
    J = J - 0.5*(wN'*WN*wN + pN'*PN*pN + wM'*WM*wM + pM'*PM*pM);
    
end

end

