%% Generate Figure 10 & 11
% This code runs for 18-24 hours.
% Reduce nExp and/or setSize for a quicker result.
clc;
clear all;
close all;

%% Experiment Parameters and Results
nExp = 100;
setSize = [10:10:90, 100:100:1000];


%% Results Varaibles
nSetSize = length(setSize);

% variables to store estimation errors
distX_geometric_SO3 = zeros(nSetSize, nExp);
distY_geometric_SO3 = zeros(nSetSize, nExp);
distX_geometric_trans = zeros(nSetSize, nExp);
distY_geometric_trans = zeros(nSetSize, nExp);

distX_prob_SO3 = zeros(nSetSize, nExp);
distY_prob_SO3 = zeros(nSetSize, nExp);
distX_prob_trans = zeros(nSetSize, nExp);
distY_prob_trans = zeros(nSetSize, nExp);

%% Distance Minimization Parameters
% Parameter Setting
alpha = 2.0;                    % translation weight
param = defaultParam;           % get default solver parameters for distance min.
% param.globalOptMethod = 2;      % activate stochastic global optimization

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark:
% The translation weight is 2.0 because the distance function in the distance
% minimization algorithm is norm(R1-R2, 'frob'), which is bounded equivalent
% to twice of geodesic distance between R1 and R2. Rotation and translation
% are equally weighted when alpha=2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% True X and Y
X_true = randSE3();  % ground truth of X
Y_true = randSE3();  % ground truth of Y


%% Noise parameters 
noiseLevel_SO3 = 0.03;              % rotation noise level in radian (std of the magnitude of angular displacement noise)
noiseLevel_trans = 0.03;            % translation noise level in user's length unit (std of translation noise)

noiseAPosition = 'left';
noiseBPosition = 'right'; 

noiseConf = 1;

%% Generate Synthetic Data
n = 1000;                            % num of measurement pairs (A,B)
[A_total, B_total] = generateABData_SE3(X_true, Y_true, n, noiseLevel_SO3, 1, noiseLevel_trans, 0.0, 'G', noiseAPosition, noiseBPosition);   % last M pairs of (A,B) are outliers.


%% Run Simulations
tic;
for j = 1:nExp
    for i = 1:length(setSize)
        %% Generate Synthetic Data
        % parameters
        N = setSize(i);                            % num of date pairs (A,B)

        % random subset of data
        idx = randperm(n);
        idx = idx(1:setSize(i));
        A = A_total(:,:,idx);
        B = B_total(:,:,idx);

        %% Solve with distance-minimization algorithm
        % Solve AX = YB with geometric stochastic global optimization
        [X_geometric1,Y_geometric1] = solveAXYB_SE3(A,B,alpha,param);
        X_geometric = X_geometric1;
        Y_geometric = Y_geometric1;


        %% Solve with probabilistic algorithm
        invSig_wN = eye(3);
        invSig_wN = repmat(invSig_wN, [1,1,n]);

        invSig_pN = eye(3);
        invSig_pN = repmat(invSig_pN, [1,1,n]);

        invSig_wM = eye(3);
        invSig_wM = repmat(invSig_wM, [1,1,n]);

        invSig_pM = eye(3);
        invSig_pM = repmat(invSig_pM, [1,1,n]);

        step_R = 0.001;
        step_p = 0.001;

        [X_prob, Y_prob, ~, ~, ~, ~] = solveAXYB_prob(A, B, X_geometric, Y_geometric, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p);


        %% Store Results
        distX_geometric_SO3(i,j) = norm(LogSO3(X_geometric(1:3,1:3) * X_true(1:3,1:3)'));
        distY_geometric_SO3(i,j) = norm(LogSO3(Y_geometric(1:3,1:3) * Y_true(1:3,1:3)'));
        distX_geometric_trans(i,j) = norm(X_geometric(1:3,4) - X_true(1:3,4));
        distY_geometric_trans(i,j) = norm(Y_geometric(1:3,4) - Y_true(1:3,4));

        distX_prob_SO3(i,j) = norm(LogSO3(X_prob(1:3,1:3) * X_true(1:3,1:3)'));
        distY_prob_SO3(i,j) = norm(LogSO3(Y_prob(1:3,1:3) * Y_true(1:3,1:3)'));
        distX_prob_trans(i,j) = norm(X_prob(1:3,4) - X_true(1:3,4));
        distY_prob_trans(i,j) = norm(Y_prob(1:3,4) - Y_true(1:3,4));

        %% Display Progress
        disp(['exp. num. = ', num2str(j), '/', num2str(nExp), ', setsize num. = ', num2str(i), '/', num2str(length(setSize))])
    end
    t = toc;
    t_left = (t/j) * (nExp-j);
    disp(['exp. num. = ', num2str(j), '/', num2str(nExp), ', time taken. = ', num2str(t), 'sec, time left = ', num2str(t_left), 'sec'])
end


%% Plot
close all

figure
subplot(2,1,1)
plot(setSize, 180/pi*mean(distX_geometric_SO3,2), setSize, 180/pi*mean(distX_prob_SO3,2));
grid on
legend('Distance minimization', 'Proposed method')
xlabel('Number of measurements')
ylabel('Rotational error (^o)')
title('Errors in rotation of X')

subplot(2,1,2)
plot(setSize, mean(distX_geometric_trans,2), setSize, mean(distX_prob_trans,2));
grid on
legend('Distance minimization', 'Proposed method')
xlabel('Number of measurements')
ylabel('Translational error')
title('Errors in translation of X')

figure
subplot(2,1,1)
plot(setSize, 180/pi*mean(distY_geometric_SO3,2), setSize, 180/pi*mean(distY_prob_SO3,2));
grid on
legend('Distance minimization', 'Proposed method')
xlabel('Number of measurements')
ylabel('Rotational error (^o)')
title('Errors in rotation of Y')

subplot(2,1,2)
plot(setSize, mean(distY_geometric_trans,2), setSize, mean(distY_prob_trans,2));
grid on
legend('Distance minimization', 'Proposed method')
xlabel('Number of measurements')
ylabel('Translational error')
title('Errors in translation of Y')


