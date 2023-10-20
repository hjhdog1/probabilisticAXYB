%% Generate Figure 6
% This code runs for 2-4 hours.
clc;
clear all;
close all;

%% Noiseless data
% number of measurement pairs
n = 20;                  

% true values of X,Y
X_true = randSE3();  % ground truth of X
Y_true = randSE3();  % ground truth of Y

% generate data
noiseAPosition = 'none';
noiseBPosition = 'none';

% generate data
[A_noiseless,B_noiseless,~] = generateABData_SE3(X_true, Y_true, n, 0, 1, 0, 0, 'G', noiseAPosition, noiseBPosition);   % last M pairs of (A,B) are outliers.


%% Noise Parameters
% parameters
noiseConf = 1;
% noiseConf = 2;
    
noiseLevel_SO3 = 0.05;         % rotation noise level in radian (std of the magnitude of angular displacement noise)
noiseLevel_trans = 0.05;            % translation noise level in user's length unit (std of translation noise)

invSig_wN = noiseLevel_SO3^-2 * eye(3);
invSig_wN = repmat(invSig_wN, [1,1,n]);

invSig_pN = noiseLevel_trans^-2 * eye(3);
invSig_pN = repmat(invSig_pN, [1,1,n]);

invSig_wM = noiseLevel_SO3^-2 * eye(3);
invSig_wM = repmat(invSig_wM, [1,1,n]);

invSig_pM = noiseLevel_trans^-2 * eye(3);
invSig_pM = repmat(invSig_pM, [1,1,n]);


%% Experiment Parameters and Results
nExp = 3000;


% variables to store solutions and errors
X = cell(nExp, 1);
Y = cell(nExp, 1);
distX_SO3 = zeros(nExp, 1);
distY_SO3 = zeros(nExp, 1);
distX_trans = zeros(nExp, 1);
distY_trans = zeros(nExp, 1);

dM = zeros(6, n, nExp);

%% Distance Minimization Parameter
% Parameter Setting
param = defaultParam;   % get default solver parameters. see instruction for more detail
% param.globalOptMethod = 2;      % stochastic global optimization with geometric local search

%% Run Experiments
tic;
for i = 1:nExp
    %% Generate Synthetic Data

    % generate data
    noiseBPosition = 'right'; 
    if noiseConf == 1
        noiseAPosition = 'left';
    elseif noiseConf == 2
        noiseAPosition = 'right';
    elseif noiseConf == 3
        noiseAPosition = 'none';
    end
    
    % generate data
    [A, B] = addNoiseSE3Data(A_noiseless, B_noiseless, noiseLevel_SO3, noiseLevel_trans, noiseAPosition, noiseBPosition, 'G');

    for j = 1:n
        M = invSE3(B_noiseless(:,:,j)) * B(:,:,j);
        dM(:,j,i) = [LogSO3(M(1:3,1:3)); M(1:3,4)];
    end
    
    %% Solve with distance-minimization algorithm for initial guess
    % Solve AX = YB with geometric stochastic global optimization
    weight = 2.0;
    [X0,Y0] = solveAXYB_SE3(A, B, weight, param);

    %% Solve with probabilistic algorithm
    if noiseConf == 3
        % probabilistic algorithm is equivalent to distance minimzation in this case
        X{i} = X0;
        Y{i} = Y0;
    else
        step_R = 1e-4;
        step_p = 1e-4;
        
        [X{i}, Y{i}, C] = solveAXYB_prob(A, B, X0, Y0, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p);
    end
    
    
    %% Display Result
    distX_SO3(i) = norm(LogSO3(X{i}(1:3,1:3) * X_true(1:3,1:3)'));
    distY_SO3(i) = norm(LogSO3(Y{i}(1:3,1:3) * Y_true(1:3,1:3)'));
    distX_trans(i) = norm(X{i}(1:3,4) - X_true(1:3,4));
    distY_trans(i) = norm(Y{i}(1:3,4) - Y_true(1:3,4));
   
    t = toc;
    disp(['====== ', num2str(i), '-th exp is over. Time taken = ', num2str(t), ', Time left = ', num2str((nExp-i) * t/i)]);
    
end

%% Results
disp(['======= Mean of errors =======']);
errMean = [mean(distX_SO3) * 180/pi;
mean(distY_SO3) * 180/pi;
mean(distX_trans);
mean(distY_trans)]

disp(['======= Std of errors =======']);
errStd = [std(distX_SO3) * 180/pi;
std(distY_SO3) * 180/pi;
std(distX_trans);
std(distY_trans)]

disp(['======= Max of errors =======']);
errMax = [max(distX_SO3) * 180/pi;
max(distY_SO3) * 180/pi;
max(distX_trans);
max(distY_trans)]


%% Compute error vectors
dX = zeros(6,nExp);
dY = zeros(6,nExp);

for i = 1:nExp
    % errors
    errX = invSE3(X_true) * X{i};
    dX(:,i) = [LogSO3(errX(1:3,1:3)); errX(1:3,4)];

    errY = invSE3(Y_true) * Y{i};
    dY(:,i) = [LogSO3(errY(1:3,1:3)); errY(1:3,4)];
end



%% Plot Analytic & Numeical Covariance

% numerical covariance
covX_n = cov(dX');
covY_n = cov(dY');

% analytic covariance
if noiseConf == 3
    [covX_a,covY_a,~] = computeUncertainty_noiseConf3(X_true, Y_true, B_noiseless, invSig_wM, invSig_pM);
    
    % estimation of analytic covariance
    [covX_est,covY_est,~] = computeUncertainty_noiseConf3(X{end}, Y{end}, B, invSig_wM, invSig_pM);
else
    C_true = zeros(4,4,n);
    for i = 1:n
        C_true(:,:,i) = A_noiseless(:,:,i) * X_true;
    end
    
    [covX_a,covY_a,~] = computeUncertainty(X_true, Y_true, C_true, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf);
    
    % estimation of analytic covariance
    [covX_est,covY_est,~] = computeUncertainty(X{end}, Y{end}, C, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf);
end

% covX_n = covY_n;
% covX_a = covY_a;
% covX_est = covY_est;


pairs = [2,3; 4,5; 2,6];
for i = 1:size(pairs,1)
    dim = pairs(i,:);

    h = figure;
    hold on
    
    scat = scatter(dX(dim(1),:), dX(dim(2),:),3,'MarkerFaceColor', 0.2*[1,1,1],'MarkerEdgeColor','none'); 
    scat.MarkerFaceAlpha = .5;
    
    covMat = covX_n(dim,dim);
    [U,S,~] = svd(covMat);
    h_ellips(1) = plotEllips([0;0], U, sqrt([S(1,1), S(2,2)]), 4);
    set(h_ellips(1), 'Color', [1,0,0])
    
    covMat = covX_a(dim,dim);
    [U,S,~] = svd(covMat);
    h_ellips(2) = plotEllips([0;0], U, sqrt([S(1,1), S(2,2)]), 3);
    set(h_ellips(2), 'Color', [0.1,0.9,0.1])

    covMat = covX_est(dim,dim);
    [U,S,~] = svd(covMat);
    h_ellips(3) = plotEllips([0;0], U, sqrt([S(1,1), S(2,2)]), 2);
    set(h_ellips(3), 'Color', [0,0,1])
    
    
	set(h,'Position',[100,100,300,350])
    set(gca,'LooseInset',get(gca,'TightInset'))
    axis equal
    
    legend(h_ellips, 'Numerical', ['Analytic (using true', newline, 'transformations)'], 'Analytic (using estimations)', 'FontSize', 12)
    
end









