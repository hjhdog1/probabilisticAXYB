function [t f_dot divide] = computeStepsize(lambda,P,Q,P0,Q0,X,Y,d,isExactStepsize);
    

    WX = skew(d(1:3));
    WY = skew(d(4:6));

    t = 0;
    f_dot = 0;
    divide = 0;
    
    if isExactStepsize
        %%%%% TODO : exact stepsize %%%%%
        t = 0;
    else
        f_dot = 0;
        divide = 0;
        for i = 1:18
            f_dot = f_dot + lambda(i) * trace(P(:,:,i)*X + Q(:,:,i)*Y) * trace(P(:,:,i)*X*WX + Q(:,:,i)*Y*WY);
%             divide = divide + 3 * abs(lambda(i)) * ( (norm(P(:,:,i)*X*WX,'fro') + norm(Q(:,:,i)*Y*WY,'fro'))^2 + (norm(P(:,:,i),'fro') + norm(Q(:,:,i),'fro')) * (norm(P(:,:,i)*X*WX^2,'fro') + norm(Q(:,:,i)*Y*WY^2,'fro')) );

        end
        
        f_dot = f_dot + trace(P0*X*WX + Q0*Y*WY);

        divide = divide + max(lambda)*(norm(WX,'fro')^2 + norm(WY,'fro')^2) + max(abs(lambda)) * sqrt(6) * sqrt(norm(WX^2,'fro')^2 + norm(WY^2,'fro')^2);
        divide = divide + sqrt(3)*(norm(P0*X*WX^2,'fro') + norm(Q0*Y*WY^2,'fro'));
        
        t = -f_dot/divide;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        step_eps = 0.01;
        if abs(t*norm(d)) < step_eps
            
            f_2dot = trace(P0*X*WX^2 + Q0*Y*WY^2);
            for i = 1:18
                f_2dot = f_2dot + lambda(i) * ( (trace(P(:,:,i)*X*WX + Q(:,:,i)*Y*WY))^2 + trace(P(:,:,i)*X + Q(:,:,i)*Y) * trace(P(:,:,i)*X*WX^2 + Q(:,:,i)*Y*WY^2) );
            end
            
            t_prev = t;
            t = -f_dot/abs(f_2dot);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            t = t * 0.5;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end

end

