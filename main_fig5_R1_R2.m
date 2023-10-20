%% Generate Figure 5
% This code runs for less than a minute (or about 10 mins when SGO is activated for distance minimization)
clc;
clear all;
close all;

%% Experiment Parameters and Results
nExp = 10;

% variables to store solutions and errors
X_geometric = cell(nExp,4);
Y_geometric = cell(nExp,4);
distX_geometric_SO3 = zeros(nExp, 4);
distY_geometric_SO3 = zeros(nExp, 4);
distX_geometric_trans = zeros(nExp, 4);
distY_geometric_trans = zeros(nExp, 4);


% X_Li = cell(nExp,4);
% Y_Li = cell(nExp,4);
% distX_Li_SO3 = zeros(nExp, 4);
% distY_Li_SO3 = zeros(nExp, 4);
% distX_Li_trans = zeros(nExp, 4);
% distY_Li_trans = zeros(nExp, 4);


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

%% Run Experiments
tic
for i = 1:nExp
    %% Generate Synthetic Data
    % parameters
    N = 20;                             % num of measurement pairs (A,B)
    noiseLevel_SO3 = 0.05;              % rotation noise level in radian (std of the magnitude of angular displacement noise)
    noiseLevel_trans = 0.05;            % translation noise level in user's length unit (std of translation noise)
    
    % true values of X,Y
    X_true = randSE3();  % ground truth of X
    Y_true = randSE3();  % ground truth of Y

    % generate data
    noiseAPosition = 'none';
    noiseBPosition = 'right';
    
    % generate data
    [A,B,M] = generateABData_SE3(X_true, Y_true, N, noiseLevel_SO3, 1, noiseLevel_trans, 0.0, 'G', noiseAPosition, noiseBPosition);   % last M pairs of (A,B) are outliers.

    
    %% Solve with distance-minimization algorithm - conf. 1: AX=YB
    % Solve AX = YB with geometric stochastic global optimization
    [X_geometric{i,1},Y_geometric{i,1}] = solveAXYB_SE3(A,B,alpha,param);
    
%     % Solve Li's method for same coordination
%     [X_Li{i,1},Y_Li{i,1}] = method.li(A,B);

    
    %% Solve with distance-minimization algorithm - conf. 2: B^-1 Y^-1 = X^-1 A^-1
    % Solve B^-1 Y^-1 = X^-1 A^-1 with geometric stochastic global optimization
    invA = invertData(A);
    invB = invertData(B);
    [invY_geometric,invX_geometric] = solveAXYB_SE3(invB,invA,alpha,param);

    X_geometric{i,2} = invSE3(invX_geometric);
    Y_geometric{i,2} = invSE3(invY_geometric);
    
%     % Solve Li's method for same coordination
%     [invY_Li, invX_Li] = method.li(invB,invA);
%     
%     X_Li{i,2} = invSE3(invX_Li);
%     Y_Li{i,2} = invSE3(invY_Li);
    
    %% Solve with distance-minimization algorithm - conf. 3: B X^-1 = Y^-1 A
    % Solve B X^-1 = Y^-1 A with geometric stochastic global optimization
    [invX_geometric,invY_geometric] = solveAXYB_SE3(B,A,alpha,param);

    X_geometric{i,3} = invSE3(invX_geometric);
    Y_geometric{i,3} = invSE3(invY_geometric);    
    
%     % Solve Li's method for same coordination
%     [invX_Li, invY_Li] = method.li(B,A);  
%     
%     X_Li{i,3} = invSE3(invX_Li);
%     Y_Li{i,3} = invSE3(invY_Li);

        
    %% Solve with distance-minimization algorithm - conf. 4: A^-1 Y = X B^-1
    % Solve A^-1 Y = X B^-1 with geometric stochastic global optimization
    invA = invertData(A);
    invB = invertData(B);
    [Y_geometric{i,4},X_geometric{i,4}] = solveAXYB_SE3(invA,invB,alpha,param);
    
%     % Solve Li's method for same coordination
%     [Y_Li{i,4}, X_Li{i,4}] = method.li(invA,invB);

    
    %% Display Result
    for j = 1:4
        distX_geometric_SO3(i,j) = norm(so3(X_geometric{i,j}(1:3,1:3) * X_true(1:3,1:3)'));
        distY_geometric_SO3(i,j) = norm(so3(Y_geometric{i,j}(1:3,1:3) * Y_true(1:3,1:3)'));
        distX_geometric_trans(i,j) = norm(X_geometric{i,j}(1:3,4) - X_true(1:3,4));
        distY_geometric_trans(i,j) = norm(Y_geometric{i,j}(1:3,4) - Y_true(1:3,4));
        
%         distX_Li_SO3(i,j) = norm(so3(X_Li{i,j}(1:3,1:3) * X_true(1:3,1:3)'));
%         distY_Li_SO3(i,j) = norm(so3(Y_Li{i,j}(1:3,1:3) * Y_true(1:3,1:3)'));
%         distX_Li_trans(i,j) = norm(X_Li{i,j}(1:3,4) - X_true(1:3,4));
%         distY_Li_trans(i,j) = norm(Y_Li{i,j}(1:3,4) - Y_true(1:3,4));
    end
    
    disp(['======= ', num2str(i), '-th experiment is done =======']);
    toc
    
end

%% Results

disp(['=========== GEOMETRIC ==========']);
disp(['======= Mean of errors =======']);
errMean = [mean(distX_geometric_SO3) * 180/pi;
mean(distY_geometric_SO3) * 180/pi;
mean(distX_geometric_trans);
mean(distY_geometric_trans)]

disp(['======= Std of errors =======']);
errStd = [std(distX_geometric_SO3) * 180/pi;
std(distY_geometric_SO3) * 180/pi;
std(distX_geometric_trans);
std(distY_geometric_trans)]

disp(['======= Max of errors =======']);
errMax = [max(distX_geometric_SO3) * 180/pi;
max(distY_geometric_SO3) * 180/pi;
max(distX_geometric_trans);
max(distY_geometric_trans)]


% disp(['================================']);
% disp(['============== Li ==============']);
% errMean_Li = [mean(distX_Li_SO3) * 180/pi;
% mean(distY_Li_SO3) * 180/pi;
% mean(distX_Li_trans);
% mean(distY_Li_trans)]
% 
% disp(['======= Std of errors =======']);
% errStd_Li = [std(distX_Li_SO3) * 180/pi;
% std(distY_Li_SO3) * 180/pi;
% std(distX_Li_trans);
% std(distY_Li_trans)]
% 
% disp(['======= Max of errors =======']);
% errMax_Li = [max(distX_Li_SO3) * 180/pi;
% max(distY_Li_SO3) * 180/pi;
% max(distX_Li_trans);
% max(distY_Li_trans)]

%% Plots
titles = {'Errors in rotation of X', 'Errors in rotation of Y';...
    'Errors in translation of X', 'Errors in translation of Y'};

close all
figure
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
for j = 1:2     % over X and Y
    for i = 1:2     % over rotation and translation
        x = 1:4;
        data = errMean(i+2*(j-1),:);
        errhigh = errStd(i+2*(j-1),:);
        errlow = errStd(i+2*(j-1),:);
%         data = [errMean(i+2*(j-1),:); errMean_Li(i+2*(j-1),:)];
%         errhigh = [errStd(i+2*(j-1),:); errStd_Li(i+2*(j-1),:)];
%         errlow = [errStd(i+2*(j-1),:); errStd_Li(i+2*(j-1),:)];

        nexttile    
        bar(x,data)
        hold on
        er = errorbar(x,data,errlow,errhigh);
        grid on

        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  

        title(titles{j,i})
        xlabel('Coordination')
        ax = gca();
        if j < 2
            ylabel('Rotational error (^o)')
            ax.YLim = [0, 3.5];
        else
            ylabel('Translational error')
            ax.YLim = [0, 0.18];
        end
    end
end


%% Box plots
close all

titles = {'Errors in rotation of X', 'Errors in rotation of Y';...
    'Errors in translation of X', 'Errors in translation of Y'};

data = {distX_geometric_SO3*180/pi, distY_geometric_SO3*180/pi, distX_geometric_trans, distY_geometric_trans};

figure
tiledlayout(2,2, 'Padding', 'none', 'TileSpacing', 'compact'); 
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
            ax.YLim = [0, 6];
        else
            ylabel('Translational error')
            ax.YLim = [0, 0.55];
        end


        pos = get(gca, 'Position');
        pos(1) = 0.08 + 0.51*(i-1);
        pos(2) = 0.08 + 0.51*(2-j);
        pos(3) = 0.405;
        pos(4) = 0.37;
        set(gca, 'Position', pos)
    end
    
end


