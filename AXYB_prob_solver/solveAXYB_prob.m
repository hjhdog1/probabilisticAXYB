function [X,Y,C,J,N,M] = solveAXYB_prob(A, B, X0, Y0, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p)

% noiselessData: either 'A' or 'B' or 'none'

if ~exist('noiselessData')
    noiselessData = 'none';
end


% initial values
X = X0;
Y = Y0;

n = size(A,3);
C = zeros(4,4,n);
for i = 1:n
    if noiselessData == 'A'
        t = 0;
    elseif noiselessData == 'A'
        t = 1;
    else
        t = (trace(invSig_wM(:,:,i)+invSig_pM(:,:,i))) / (trace(invSig_wN(:,:,i)+invSig_pN(:,:,i)+invSig_wM(:,:,i)+invSig_pM(:,:,i)));
    end
    C(:,:,i) = interpolateSE3(A(:,:,i)*X, Y*B(:,:,i), t);
end


% mat2cell
Ac = num2cell(A,[1,2]);
Bc = num2cell(B,[1,2]);
C = num2cell(C,[1,2]);
invSig_wNc = num2cell(invSig_wN,[1,2]);
invSig_pNc = num2cell(invSig_pN,[1,2]);
invSig_wMc = num2cell(invSig_wM,[1,2]);
invSig_pMc = num2cell(invSig_pM,[1,2]);

% iteration
for i = 1:5000
    [GwX, GwY, GwC, GqX, GqY, GqC, J] = gradLogL(X, Y, Ac, Bc, C, invSig_wNc, invSig_pNc, invSig_wMc, invSig_pMc, noiseConf);
%     [X,Y,C] = updateXYC(X, Y, C, GwX, GwY, GwC, GqX, GqY, GqC, step_R, step_p, n);
    [X,Y,C] = updateXYC2(X, Y, C, GwX, GwY, GwC, GqX, GqY, GqC, step_R, step_p, n);
%     J = computeLogL(X, Y, A, B, C, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf);
    
    if mod(i,200) == 0
%         disp(GwX(1)*1e6)
%         disp(i)
%         disp(J)
        % SO(3) projection
        X(1:3,1:3) = maxTraceFX(X(1:3,1:3)');
        Y(1:3,1:3) = maxTraceFX(Y(1:3,1:3)');
        for j = 1:n
%             C(1:3,1:3, j) = maxTraceFX(C(1:3,1:3, j)');
            C{j}(1:3,1:3) = maxTraceFX(C{j}(1:3,1:3)');
        end
    end
end

% cell2mat
C = cell2mat(C);

% N and M
N = zeros(4,4,n);
M = zeros(4,4,n);
for i = 1:n
    if noiseConf == 1
        N(:,:,i) = A(:,:,i) * X * invSE3(C(:,:,i));
    else
        N(:,:,i) = X * invSE3(C(:,:,i)) * A(:,:,i);
    end
    
    M(:,:,i) = invSE3(C(:,:,i)) * Y * B(:,:,i);
    
    
end

% final cost
J = computeLogL(X, Y, A, B, C, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf);

end

