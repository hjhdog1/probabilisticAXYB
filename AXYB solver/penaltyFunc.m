function f = penaltyFunc(x,A,B,alpha)
    N = size(A,3);
    
    [X Y] = v2XY_SE3(x);
    
    mu1 = 1;
    mu2 = alpha;
    mu3 = 1e6;
    mu4 = 1e6;
    
    f = 0;
    
    for i = 1:N
        diff = A(:,:,i)*X - Y*B(:,:,i);
        fnorm_diff = sum(sum(diff(1:3,1:3).*diff(1:3,1:3)));
        f = f + mu1*fnorm_diff;
        
        fnorm_diff = sum(diff(1:3,4).*diff(1:3,4));
        f = f + mu2*fnorm_diff;
        
        diff = X(1:3,1:3)*X(1:3,1:3)' - eye(3);
        fnorm_diff = sum(sum(diff.*diff));
        f = f + mu3*fnorm_diff;
        
        diff = Y(1:3,1:3)*Y(1:3,1:3)' - eye(3);
        fnorm_diff = sum(sum(diff.*diff));
        f = f + mu4*fnorm_diff;
        
    end

end

