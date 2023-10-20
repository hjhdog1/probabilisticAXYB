function [X,Y,J] = solveAXYB_prob_noiselessA(A, B, X0, Y0, invSig_wM, invSig_pM, step_R, step_p)
% A is noiseless and B is noisy (M = X^{-1}A^{-1}YB)

% initial values
X = X0;
Y = Y0;

n = size(A,3);


% mat2cell
Ac = num2cell(A,[1,2]);
Bc = num2cell(B,[1,2]);
invSig_wMc = num2cell(invSig_wM,[1,2]);
invSig_pMc = num2cell(invSig_pM,[1,2]);

% iteration
for i = 1:5000
    [GwX, GwY, GqX, GqY, J] = gradLogL_noiselessA(X, Y, Ac, Bc, invSig_wMc, invSig_pMc);
    [X,Y] = updateXYC(X, Y, zeros(0,0,0), GwX, GwY, zeros(0,0,0), GqX, GqY, zeros(0,0,0), step_R, step_p, n);
    
    if mod(i,200) == 0
%         disp(J)
        % SO(3) projection
        X(1:3,1:3) = maxTraceFX(X(1:3,1:3)');
        Y(1:3,1:3) = maxTraceFX(Y(1:3,1:3)');
    end
end


% final cost
J = computeLogL_noiselessA(X, Y, A, B, invSig_wM, invSig_pM);

end

