function AA = buildH11(A,B)
% A,B : SO(3)
    N = size(A,3);
    AA = zeros(18);
    
    for i = 1:N
        Q_hat = zeros(9);
        P_hat = zeros(9);
        
        for j = 1:3
            Q_hat(1+3*(j-1):3*j, 1+3*(j-1):3*j) = B(:,:,i)';
%             P_hat([0+j 3+j 6+j], 1+3*(j-1):3*j) = A(:,:,i)';
            P_hat([0+j 3+j 6+j], [0+j 3+j 6+j]) = A(:,:,i)';
        end
        AA(1:9, 10:18) = AA(1:9, 10:18) + Q_hat*P_hat;
    end
    AA = -1*AA;
    AA = AA+AA';
end

