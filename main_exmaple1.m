%% Example 1: running calibration on synthetic data
clc;
clear all;
close all;

%% Generate Synthetic Data
X_true = randSE3();
Y_true = randSE3();

n = 20;                             % num of measurement pairs (A,B)

A = zeros(4,4,n);
B = zeros(4,4,n);
for i = 1:n
    A(:,:,i) = randSE3();
    B(:,:,i) = invSE3(Y_true) * A(:,:,i) * X_true;
end


%% Set Noise configuration and noisy properties
noiseConf = 1;  % This can be either 1, 2 or 3 depending on your system noise. See the paper for details.

% Noise covariance matrices of position and translation of A and B
Sig_wN = repmat(eye(3), [1,1,n]);   
Sig_pN = repmat(eye(3), [1,1,n]);
Sig_wM = repmat(eye(3), [1,1,n]);
Sig_pM = repmat(eye(3), [1,1,n]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matching noiseConf and noise covariances with your system noise will give
% you the best calibration result. % In this particular example, however,
% the data is noiseless; thus, any noiseConf and noise covariances will
% give the same calibration result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run Calibration
% Noise covariances are inverted and fed to solveAXYB_prob() function.
invSig_wN = invert_cov(Sig_wN);
invSig_pN = invert_cov(Sig_pN);
invSig_wM = invert_cov(Sig_wM);
invSig_pM = invert_cov(Sig_pM);

% A reasonable initial guesse is needed. One can be obtained by any
% existing method. Here we assume it is available.
X0 = X_true * [LargeSO3(0.5*randn(3,1)), 0.5*randn(3,1); 0,0,0,1];
Y0 = Y_true * [LargeSO3(0.5*randn(3,1)), 0.5*randn(3,1); 0,0,0,1];

% Solving
if noiseConf == 1 || noiseConf == 2
    [X_est, Y_est] = solveAXYB_prob(A, B, X0, Y0, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, 0.05, 0.05);
elseif noiseConf == 3
    [X_est, Y_est] = solveAXYB_prob_noiselessA(A, B, X0, Y0, invSig_wM, invSig_pM, 0.05, 0.05);
end


%% Print Errors
disp(X_true - X_est)
disp(Y_true - Y_est)