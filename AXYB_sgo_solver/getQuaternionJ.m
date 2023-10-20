function [J qx qy] = getQuaternionJ(C_i,signs)
        % 초기값의 결정
        N = length(signs);
        C = zeros(4);

        for i = 1:N
            C = C - signs(i) * C_i(:,:,i);
        end
        S = [N*eye(4) C; C' N*eye(4)];
        CC = C'*C;
        [vec val] = eig(CC);
        val = diag(val);

        lambda_plus = N*ones(4,1) + sqrt(val);
        lambda_minus = N*ones(4,1) - sqrt(val);
        lambda = [lambda_plus ; lambda_minus];

        lambda(lambda < 0) = max(lambda);

        [minLambda minIdx] = min(lambda);
        isMinus = 0;
        if minIdx > 4
            minIdx = minIdx - 4;
            isMinus = 1;
        end

        qy = vec(:,minIdx)*sign(vec(1,minIdx));
        qx = 1/(minLambda-N)*C*qy;
        
        J = [qx;qy]' * S * [qx;qy];
end

