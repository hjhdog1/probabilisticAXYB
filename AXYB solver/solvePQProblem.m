function [X,Y,J] = solvePQProblem(lambda,P,Q,P0,Q0,c,X0,Y0,maxIter,ftol,isNewtonMethod,isExactStepsize)
% X,Y : SO(3)

    slowConvergence = 0;
    step_eps = 0.0001;

    X = X0;
    Y = Y0;
    J_prev = computeCost_SE3_PQ(lambda,P,Q,P0,Q0,c,X,Y);
    for i = 1:maxIter
        
        % find direction and stepsize
        d = computeDirection_PQ(lambda,P,Q,P0,Q0,X,Y,isNewtonMethod);
        t = 1;
        if ~slowConvergence
            t = computeStepsize(lambda,P,Q,P0,Q0,X,Y,d,isExactStepsize);
        end
               
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if abs(t*norm(d)) < step_eps
%             slowConvergence = 1;
%             isNewtonMethod = 1;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if slowConvergence && (abs(norm(d)) > 1)
%             slowConvergence = 0;
%             isNewtonMethod = 0;
%             continue;
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        % update
        X = X*LargeSO3(d(1:3)*t);
        Y = Y*LargeSO3(d(4:6)*t);
        
        % compute J
        J = computeCost_SE3_PQ(lambda,P,Q,P0,Q0,c,X,Y);
        
%         if t*norm(d)/J < ftol
        if t*norm(d) < ftol
            break;
        end
 
        
    end

end

