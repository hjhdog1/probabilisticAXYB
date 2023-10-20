%% Generate Figure 7 & 8
% This code runs for 10 minutes (or about 15 mins when SGO is activated for distance minimization)
% change noiseConf to choose noise configuration between 1 and 2
clc;
clear all;
close all;

%% Set Noise configuration
noiseConf = 1;  % for figure 7
% noiseConf = 2;  % for figure 8

%% Experiment Parameters and Result Variables
methods = {'Park', 'Li', 'Shah', 'Tabb_Zc2', 'geometric', 'prob'};
nMethods = length(methods);

nExp = 100;

% variables to store estimation errors
distX_SO3 = zeros(nExp, nMethods);
distY_SO3 = zeros(nExp, nMethods);
distX_trans = zeros(nExp, nMethods);
distY_trans = zeros(nExp, nMethods);

% variables to store estimations
X = cell(nExp, nMethods);
Y = cell(nExp, nMethods);


%% Geometric Distance Minimization Parameters
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



%% Run Calibrations
timer = tic;
for i = 1:nExp
    %% Generate Synthetic Data
    % parameters
    nMeas = 20;                             % num of measurement pairs (A,B)
    noiseLevel_SO3 = 0.05;              % rotation noise level in radian (std of the magnitude of angular displacement noise)
    noiseLevel_trans = 0.05;            % translation noise level in user's length unit (std of translation noise)
    
    % true values of X,Y
    X_true{i} = randSE3();  % ground truth of X
    Y_true{i} = randSE3();  % ground truth of Y

    % generate data
    noiseBPosition = 'right'; 
    if noiseConf == 1
        noiseAPosition = 'left';
    elseif noiseConf == 2
        noiseAPosition = 'right';
    end
    
    
    % generate data
    [A,B] = generateABData_SE3(X_true{i}, Y_true{i}, nMeas, noiseLevel_SO3, 1, noiseLevel_trans, 0.0, 'G', noiseAPosition, noiseBPosition);   % last M pairs of (A,B) are outliers.
    [A,B] = randomSorting(A,B); % random-sorting data to mix outliers
    
    % generate inverted data
    invA = invertData(A);
    invB = invertData(B);

    % run methods
    for j = 1:nMethods
        switch methods{j}
            case 'Park'
                X{i,j} = method.navy_calibration(A, invB);
                Y{i,j} = method.navy_calibration(invA, B);
            case 'Li'                
                [X{i,j}, Y{i,j}] = method.li(A,B);
            case 'Tabb_Zc2'
%                 [~,handEyeHT.TabbZc1,base2grid.TabbZc2,handEyeHT.TabbZc2,time.TabbZc1,time.TabbZc2] = method.AXZB(BasePoseIntcpCoords,gridPoseInCameraCoords,initial,Config) ;
                Config=H2EConfiguration;
                Config.SolAlgo='levenberg_marquardt';
                initial=[ rotm2quat(X{i,j-1}(1:3,1:3)) X{i,j-1}(1:3,4)'  rotm2quat(Y{i,j-1}(1:3,1:3)), Y{i,j-1}(1:3,4)'];
%                 initial=[ rot2quat(X{i,j-1}(1:3,1:3)) X{i,j-1}(1:3,4)  rot2quat(Y{i,j-1}(1:3,1:3)), Y{i,j-1}(1:3,4)];
                [~,~,X{i,j}, Y{i,j},~,~] = method.AXZB(A,B,initial,Config) ;                
            case 'Shah'
                [X{i,j}, Y{i,j}] = method.HandeyeShah(A,B) ; %pose of grid relative to the camera (Hcam2grid2)
            case 'geometric'                
                [X{i,j}, Y{i,j}] = solveAXYB_SE3(A,B,alpha,param);
            case 'prob'    
                Eyes = repmat(eye(3), [1,1,nMeas]);
                invSig_wN = Eyes;
                invSig_pN = Eyes;
                invSig_wM = Eyes;
                invSig_pM = Eyes;

                step_R = 0.01;
                step_p = 0.01;

                [X{i,j}, Y{i,j}, C, J, N, M] = solveAXYB_prob(A, B, X{i,j-1}, Y{i,j-1}, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p);
        end
        
    end

    
    %% Store Results
    for j = 1:nMethods
        distX_SO3(i,j) = norm(LogSO3(X{i,j}(1:3,1:3) * X_true{i}(1:3,1:3)'));
        distY_SO3(i,j) = norm(LogSO3(Y{i,j}(1:3,1:3) * Y_true{i}(1:3,1:3)'));
        distX_trans(i,j) = norm(X{i,j}(1:3,4) - X_true{i}(1:3,4));
        distY_trans(i,j) = norm(Y{i,j}(1:3,4) - Y_true{i}(1:3,4));
    end

    %% Display Time
    t = toc(timer);
    t_left = (t/i) * (nExp-i);
    disp(['exp. num. = ', num2str(i), '/', num2str(nExp), ', time taken. = ', num2str(t), 'sec, time left = ', num2str(t_left), 'sec'])


end

errMean_X_SO3 = mean(distX_SO3) * 180/pi;
errMean_Y_SO3 = mean(distY_SO3) * 180/pi;
errMean_X_trans = mean(distX_trans);
errMean_Y_trans = mean(distY_trans);

errMean = [errMean_X_SO3; errMean_Y_SO3; errMean_X_trans; errMean_Y_trans];

errStd_X_SO3 = std(distX_SO3) * 180/pi;
errStd_Y_SO3 = std(distY_SO3) * 180/pi;
errStd_X_trans = std(distX_trans);
errStd_Y_trans = std(distY_trans);

errStd = [errStd_X_SO3; errStd_Y_SO3; errStd_X_trans; errStd_Y_trans];


%% Plot Results - Box plots
titles = {'Errors in rotation of X', 'Errors in rotation of Y';...
    'Errors in translation of X', 'Errors in translation of Y'};

data = {distX_SO3*180/pi, distY_SO3*180/pi, distX_trans, distY_trans};

figure('Position', [100, 100, 700, 500])
% figure
for j = 1:2     % over X and Y
    for i = 1:2     % over rotation and translation
        x = 1:5;

        subplot(2,2,2*(j-1)+i)

        boxplot(data{2*(j-1)+i},'Whisker',100)
        grid on

        title(titles{j,i})
        ax = gca;
        if j < 2
            ylabel('Rotational error (^o)')
            ax.YLim = [0, 6.2];
        else
            ylabel('Translational error')
            ax.YLim = [0, 0.45];
        end

%         xticklabels({'AX=XB', 'Li', 'Shah', 'Dist. Min.', 'Proposed'})
        xticklabels({'AX=XB', 'Li', 'Shah', 'Tabb', 'Dist. Min.', 'Proposed'})
        xtickangle(30)
        ax = gca;
        xrule = ax.XAxis;
        xrule.FontSize = 8.5;
        ax.XRuler.TickLabelGapOffset = -2;

        pos = get(gca, 'Position');
        pos(1) = 0.08 + 0.51*(i-1);
        pos(2) = 0.08 + 0.51*(2-j);
        pos(3) = 0.405;
        pos(4) = 0.37;
        set(gca, 'Position', pos)
    end
    
end
