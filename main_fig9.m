%% Generate Figure 9
% This code runs for 5-10 minutes (or about 8-15 mins when SGO is activated for distance minimization)
clc;
clear all;
close all;


%% Experiment Parameters and Results
nExp = 100;

% variables to store estimation errors
distX_geometric_SO3 = zeros(1, nExp);
distY_geometric_SO3 = zeros(1, nExp);
distX_geometric_trans = zeros(1, nExp);
distY_geometric_trans = zeros(1, nExp);

distX_prob_SO3 = zeros(1, nExp);
distY_prob_SO3 = zeros(1, nExp);
distX_prob_trans = zeros(1, nExp);
distY_prob_trans = zeros(1, nExp);

diffX_SO3 = zeros(1, nExp);
diffY_SO3 = zeros(1, nExp);
diffX_trans = zeros(1, nExp);
diffY_trans = zeros(1, nExp);

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


%% Run Calibrations
tic;
for i = 1:nExp
    %% Generate Synthetic Data
    % parameters
    n = 20;                            % num of date pairs (A,B)
    noiseLevel_SO3 = 0.05;         % rotation noise level in radian (std of the magnitude of angular displacement noise)
    noiseLevel_trans = 0.05;            % translation noise level in user's length unit (std of translation noise)
    
    % true values of X,Y
    X_true = randSE3();  % ground truth of X
    Y_true = randSE3();  % ground truth of Y

    % generate data
    noiseBPosition = 'right'; 
    noiseAPosition = 'none';
    
    
    % generate data
    [A,B] = generateABData_SE3(X_true, Y_true, n, noiseLevel_SO3, 1, noiseLevel_trans, 0.0, 'G', noiseAPosition, noiseBPosition);   % last M pairs of (A,B) are outliers.
    [A,B] = randomSorting(A,B); % random-sorting data to mix outliers


    %% Solve with distance-minimization algorithm
    % Solve AX = YB with geometric stochastic global optimization
    [X_geometric1,Y_geometric1] = solveAXYB_SE3(A,B,alpha,param);
    X_geometric = X_geometric1;
    Y_geometric = Y_geometric1;

    
    %% Solve with probabilistic algorithm
    invSig_wM = eye(3);
    invSig_wM = repmat(invSig_wM, [1,1,n]);

    invSig_pM = eye(3);
    invSig_pM = repmat(invSig_pM, [1,1,n]);

    
    X_prob = addNoiseSE3(X_geometric, 0.02, 0.02, 'right', 'G');
    Y_prob = addNoiseSE3(Y_geometric, 0.02, 0.02, 'right', 'G');

    step_R = 0.01;
    step_p = 0.01;
    [X_prob, Y_prob, ~] = solveAXYB_prob_noiselessA(A, B, X_prob, Y_prob, invSig_wM, invSig_pM, step_R, step_p);

    step_R = 0.0001;
    step_p = 0.0001;
    [X_prob, Y_prob, J] = solveAXYB_prob_noiselessA(A, B, X_prob, Y_prob, invSig_wM, invSig_pM, step_R, step_p);

    
    %% Store Results
    distX_geometric_SO3(i) = norm(LogSO3(X_geometric(1:3,1:3) * X_true(1:3,1:3)'));
    distY_geometric_SO3(i) = norm(LogSO3(Y_geometric(1:3,1:3) * Y_true(1:3,1:3)'));
    distX_geometric_trans(i) = norm(X_geometric(1:3,4) - X_true(1:3,4));
    distY_geometric_trans(i) = norm(Y_geometric(1:3,4) - Y_true(1:3,4));

    distX_prob_SO3(i) = norm(LogSO3(X_prob(1:3,1:3) * X_true(1:3,1:3)'));
    distY_prob_SO3(i) = norm(LogSO3(Y_prob(1:3,1:3) * Y_true(1:3,1:3)'));
    distX_prob_trans(i) = norm(X_prob(1:3,4) - X_true(1:3,4));
    distY_prob_trans(i) = norm(Y_prob(1:3,4) - Y_true(1:3,4));
    
    diffX_SO3(i) = norm(LogSO3(X_geometric(1:3,1:3) * X_prob(1:3,1:3)'));
    diffY_SO3(i) = norm(LogSO3(Y_geometric(1:3,1:3) * Y_prob(1:3,1:3)'));
    diffX_trans(i) = norm(X_geometric(1:3,4) - X_prob(1:3,4));
    diffY_trans(i) = norm(Y_geometric(1:3,4) - Y_prob(1:3,4));

    %% Display Time
    t = toc;
    t_left = (t/i) * (nExp-i);
    disp(['exp. num. = ', num2str(i), '/', num2str(nExp), ', time taken. = ', num2str(t), 'sec, time left = ', num2str(t_left), 'sec'])


end

errMean = [mean(distX_geometric_SO3) * 180/pi, mean(distX_prob_SO3) * 180/pi;
mean(distY_geometric_SO3) * 180/pi, mean(distY_prob_SO3) * 180/pi;
mean(distX_geometric_trans), mean(distX_prob_trans);
mean(distY_geometric_trans), mean(distY_prob_trans)]

errStd = [std(distX_geometric_SO3) * 180/pi, std(distX_prob_SO3) * 180/pi;
std(distY_geometric_SO3) * 180/pi, std(distY_prob_SO3) * 180/pi;
std(distX_geometric_trans), std(distX_prob_trans);
std(distY_geometric_trans), std(distY_prob_trans)]

%% Plot Results
close all

titles = {'Errors in rotation of X', 'Errors in rotation of Y';...
    'Errors in translation of X', 'Errors in translation of Y'};


figure
% tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
for j = 1:2     % over X and Y
    for i = 1:2     % over rotation and translation
        x = 1:2;
        data = errMean(i+2*(j-1),:);
        errhigh = errStd(i+2*(j-1),:);
        errlow = errStd(i+2*(j-1),:);

        subplot(2,2,2*(j-1)+i)
%         nexttile
        bar(x,data)
        hold on
        er = errorbar(x,data,errlow,errhigh);
        grid on

        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  

        title(titles{j,i})
%         xlabel('Methods')
        ax = gca();
        if j < 2
            ylabel('Rotational error (^o)')
            ax.YLim = [0, 2.0];
        else
            ylabel('Translational error')
            ax.YLim = [0, 0.1];
        end
%         xticklabels({'distance minimization', 'proposed method'})
        xticklabels({'   Distance\newlineminimization', 'Proposed\newline method'})
%         ax = gca();
%         ax.XTickLabel = tickLabels;
        pos = get(gca, 'Position');
        pos(1) = 0.08 + 0.51*(i-1);
        pos(2) = 0.08 + 0.51*(2-j);
        pos(3) = 0.405;
        pos(4) = 0.37;
        set(gca, 'Position', pos)
    end
    
end

mean([diffX_SO3', diffY_SO3', diffX_trans', diffY_trans'])
max([diffX_SO3', diffY_SO3', diffX_trans', diffY_trans']) ./ errMean(:,1)' * 100
min([diffX_SO3', diffY_SO3', diffX_trans', diffY_trans']) ./ errMean(:,1)' * 100


%% Plot Results 3 - Box plots
close all

titles = {'Errors in rotation of X', 'Errors in rotation of Y';...
    'Errors in translation of X', 'Errors in translation of Y'};

data = {[distX_geometric_SO3', distX_prob_SO3']*180/pi, [distY_geometric_SO3', distY_prob_SO3']*180/pi,...
        [distX_geometric_trans', distX_prob_trans'], [distY_geometric_trans', distY_prob_trans']};

figure('Position', [100, 100, 700, 500])
% figure
for j = 1:2     % over X and Y
    for i = 1:2     % over rotation and translation
        x = 1:2;

        subplot(2,2,2*(j-1)+i)

        bp = boxplot(data{2*(j-1)+i},'Whisker',100);
        grid on

        title(titles{j,i})
        ax = gca();
        if j < 2
            ylabel('Rotational error (^o)')
            ax.YLim = [0, 2.0];
        else
            ylabel('Translational error')
            ax.YLim = [0, 0.06];
        end

%         xticklabels({'   Distance\newlineminimization', 'Proposed\newline method'})
        xticklabels({'Dist. Min.', 'Proposed'})

        pos = get(gca, 'Position');
        pos(1) = 0.08 + 0.51*(i-1);
        pos(2) = 0.08 + 0.51*(2-j);
        pos(3) = 0.405;
        pos(4) = 0.37;
        set(gca, 'Position', pos)
    end
    
end
