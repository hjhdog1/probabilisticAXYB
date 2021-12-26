function [X,Y,dummy,methodName] = solveAXYB_SE3(A,B,alpha,param)
% param contains method, globalOptMethod, maxIter,ftol,sampleNumPerIter
% param.method : 1 - steepest, stepsize est
%                2 - newton, stepsize est
%                3 - steepest, exact stepsize
%                4 - newton, exact stepsize
% param.globalOptMethod : 0 - no global opt
%                         1 - multistart
%                         2 - gamma-clustering
%                        -1 - quaternion-based
%                        -2 - quaternion-based with genetic algorithm
% param.sampleNumPerIter : sampling number per iteration in global optimization
dummy = struct;

N = size(A,3);
isSO3Problem = false;
if size(A,1) == 3
    
    isSO3Problem = true;
    
    alpha = 0;
    
    A_SO3 = A;
    B_SO3 = B;
    
    A = zeros(4,4,N);
    B = zeros(4,4,N);
    
    for i = 1:N
        A(1:3,1:3,i) = A_SO3(:,:,i);
        A(4,4,i) = 1;
        B(1:3,1:3,i) = B_SO3(:,:,i);
        B(4,4,i) = 1;
    end    
end

isNewtonMethod = 0;
isExactStepsize = 0;

if param.method == 1
    isNewtonMethod = 0;
    isExactStepsize = 0;
end
if param.method == 2
    isNewtonMethod = 1;
    isExactStepsize = 0;
end
if param.method == 3
    isNewtonMethod = 0;
    isExactStepsize = 1;
end
if param.method == 4
    isNewtonMethod = 1;
    isExactStepsize = 1;
end


methodName = 'geometric optimization';

X1 = randSO3;
Y1 = randSO3;

%% Quaternion based method
if param.globalOptMethod == -1
    [H,f,r] = getHfr(A,B,alpha);


    Q_A = zeros(4,4,N);
    W_B = zeros(4,4,N);
    C_i = zeros(4,4,N);
    for i = 1:N
        Q_A(:,:,i) = q2Q(SO3toQuart(A(1:3,1:3,i)));
        W_B(:,:,i) = q2W(SO3toQuart(B(1:3,1:3,i)));
        C_i(:,:,i) = Q_A(:,:,i)' * W_B(:,:,i);
    end

    %%%%%%%%%%%%%%%%%% Book Keeping %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1 = quaternionBookKepping(A) .* quaternionBookKepping(B);
    x2 = -x1;
    
    [J1, qx1, qy1] = getQuaternionJ(C_i, x1);
    [J2, qx2, qy2] = getQuaternionJ(C_i, x2);
    
    if J1 < J2
        qx = qx1;
        qy = qy1;
    else
        qx = qx2;
        qy = qy2;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = Quart2SO3(qx);
    Y = Quart2SO3(qy);
    [X,Y] = getOptSE3(X,Y,H,f);
    
    % 최적화
    options = optimoptions('fminunc','MaxIter',param.maxIter,'TolFun',param.ftol,'TolX',0,'MaxFunEvals',inf);
    
    x0 = XY2v_SE3(X,Y);
    x = fminunc(@(x)penaltyFunc(x,A,B,alpha), x0, options);
    [X,Y] = v2XY_SE3(x);

    X(1:3,1:3) = maxTraceFX(X(1:3,1:3)');
    Y(1:3,1:3) = maxTraceFX(Y(1:3,1:3)');
    
    methodName = 'local quaternion (with closed form initial guess)';

    J = computeCost_SE3_reduced(H,f,r,X,Y);
    dummy.XY_minima.X = X(1:3,1:3);
    dummy.XY_minima.Y = Y(1:3,1:3);
    dummy.XY_minima.J = J;
    
end

%% Quaternion based method-GA
if param.globalOptMethod == -2
    [H,f,r] = getHfr(A,B,alpha);


    Q_A = zeros(4,4,N);
    W_B = zeros(4,4,N);
    C_i = zeros(4,4,N);
    for i = 1:N
        Q_A(:,:,i) = q2Q(SO3toQuart(A(1:3,1:3,i)));
        W_B(:,:,i) = q2W(SO3toQuart(B(1:3,1:3,i)));
        C_i(:,:,i) = Q_A(:,:,i)' * W_B(:,:,i);
    end
  
    %%%%%%%%%%%%%%%%% Genetic Algorithm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fun = @(x)getQuaternionJ(C_i, (-ones(1,N)).^x );
    x = ga(fun,N,[],[],[],[],-0.5*ones(N,1),1.5*ones(N,1),[],1:N);
    [~, qx, qy] = fun(x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X = Quart2SO3(qx);
    Y = Quart2SO3(qy);
    [X,Y] = getOptSE3(X,Y,H,f);
    
    % 최적화
    options = optimoptions('fminunc','MaxIter',param.maxIter,'TolFun',param.ftol,'TolX',0,'MaxFunEvals',inf);
    
    x0 = XY2v_SE3(X,Y);
    x = fminunc(@(x)penaltyFunc(x,A,B,alpha), x0, options);
    [X,Y] = v2XY_SE3(x);

    X(1:3,1:3) = maxTraceFX(X(1:3,1:3)');
    Y(1:3,1:3) = maxTraceFX(Y(1:3,1:3)');
    
    methodName = 'local quaternion (with initial guess from genetic algorithm)';

    J = computeCost_SE3_reduced(H,f,r,X,Y);
    dummy.XY_minima.X = X(1:3,1:3);
    dummy.XY_minima.Y = Y(1:3,1:3);
    dummy.XY_minima.J = J;
    
end

%% Local Optimization with closed-from initial guess
if param.globalOptMethod == 0
    [H,f,r] = getHfr(A,B,alpha);
    [H_subOpt,f_subOpt,r_subOpt] = getHfr_subOpt(H,f,r);
    [lambda,P,Q,P0,Q0,c] = getPQ(H_subOpt,f_subOpt,r_subOpt);

    [X0, Y0] = findInitGuessUsingAXXB2(A(1:3,1:3,:),B(1:3,1:3,:));
    
    [X,Y,J] = solvePQProblem(lambda,P,Q,P0,Q0,c,X0(1:3,1:3),Y0(1:3,1:3),param.maxIter,param.ftol,isNewtonMethod,isExactStepsize);
    
    dummy.XY_minima.X = X;
    dummy.XY_minima.Y = Y;
    dummy.XY_minima.J = J;
    
    [X,Y] = getOptSE3(X,Y,H,f);
    methodName = 'local geometric (with closed form initial guess)';
    
end

%% Local Optimization with random initial guess
if param.globalOptMethod == -20
    [H,f,r] = getHfr(A,B,alpha);
    [H_subOpt,f_subOpt,r_subOpt] = getHfr_subOpt(H,f,r);
    [lambda,P,Q,P0,Q0,c] = getPQ(H_subOpt,f_subOpt,r_subOpt);
 
    [X0, Y0] = findEigenInitXY(P,Q);
    
    [X,Y,J] = solvePQProblem(lambda,P,Q,P0,Q0,c,X0(1:3,1:3),Y0(1:3,1:3),param.maxIter,param.ftol,isNewtonMethod,isExactStepsize);
    dummy.XY_minima.X = X;
    dummy.XY_minima.Y = Y;
    dummy.XY_minima.J = J;
    
    [X,Y] = getOptSE3(X,Y,H,f);
    methodName = 'global geometric (with random initial guess)';
    
end

%% Multistart
if param.globalOptMethod == 1
    [H,f,r] = getHfr(A,B,alpha);
    [H_subOpt,f_subOpt,r_subOpt] = getHfr_subOpt(H,f,r);
    [lambda,P,Q,P0,Q0,c] = getPQ(H_subOpt,f_subOpt,r_subOpt);
        
    XYsample = struct;
    XY_minima = [];
    
    k = 1;
    isStopping = 0;
    
    while ~isStopping
        XYsample(k).X = randMultipleSO3(param.sampleNumPerIter);
        XYsample(k).Y = randMultipleSO3(param.sampleNumPerIter);
        
        for i = 1:param.sampleNumPerIter
            [X,Y,J] = solvePQProblem(lambda,P,Q,P0,Q0,c,XYsample(k).X(:,:,i),XYsample(k).Y(:,:,i),param.maxIter,param.ftol,isNewtonMethod,isExactStepsize);
            XY_minima = addXY(XY_minima,X,Y,J);
        end
        
        w = length(XY_minima);
        isStopping = stopping(w,k*param.sampleNumPerIter);
            
        k = k+1;
    end
    
    [X,Y] = findBestXY(XY_minima);
    [X,Y] = getOptSE3(X,Y,H,f);
    
    dummy.XY_minima = XY_minima;
    methodName = 'global geometric (with random initial guess)';
    
end


%% Gamma-Clustering
if param.globalOptMethod == 2
    [H,f,r] = getHfr(A,B,alpha);
    [H_subOpt,f_subOpt,r_subOpt] = getHfr_subOpt(H,f,r);
    [lambda,P,Q,P0,Q0,c] = getPQ(H_subOpt,f_subOpt,r_subOpt);
        
    XYsample = struct;
    XY_minima = [];
    
    k = 1;
    isStopping = 0;
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%% Usign Initial Guess %%%%%%%%%%%%%%%%%
    [X0, Y0] = findInitGuessUsingAXXB2(A(1:3,1:3,:),B(1:3,1:3,:));
  
    [X,Y,J] = solvePQProblem(lambda,P,Q,P0,Q0,c,X0(1:3,1:3),Y0(1:3,1:3),param.maxIter,param.ftol,isNewtonMethod,isExactStepsize);
    XY_minima = addXY(XY_minima,X,Y,J);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    Js = [];
%     Loop_Iter = 0;
    while ~isStopping
%         Loop_Iter = Loop_Iter + 1;
        
        XYsample(k).X = randMultipleSO3(param.sampleNumPerIter);
        XYsample(k).Y = randMultipleSO3(param.sampleNumPerIter);
        
        Js_cur = zeros(1,param.sampleNumPerIter);
        for i = 1:param.sampleNumPerIter
            Js_cur(i) = computeCost_SE3_PQ(lambda,P,Q,P0,Q0,c,XYsample(k).X(:,:,i),XYsample(k).Y(:,:,i));
        end
        Js = [Js Js_cur];
        Js = sort(Js);
        levelcut = getPropositionalLevel(Js, param.gamma);
        
        for i = 1:param.sampleNumPerIter
            if Js_cur(i) >= levelcut
                continue;
            end
            
            [X,Y,J] = solvePQProblem(lambda,P,Q,P0,Q0,c,XYsample(k).X(:,:,i),XYsample(k).Y(:,:,i),param.maxIter,param.ftol,isNewtonMethod,isExactStepsize);
            XY_minima = addXY(XY_minima,X,Y,J);
        end
        
        w = length(XY_minima);
        isStopping = stopping(w,param.gamma*k*param.sampleNumPerIter);
            
        k = k+1;
        
%                 
%         if Loop_Iter > 1
%            break; 
%         end
    end
    
    [X,Y] = findBestXY(XY_minima);
    [X,Y] = getOptSE3(X,Y,H,f);
        
    dummy.XY_minima = XY_minima;
    methodName = 'global geometric (with random initial guess)';
    
end


%% Cutting out unnecessary translation part
if isSO3Problem
    X = X(1:3,1:3);
    Y = Y(1:3,1:3);
end



end

