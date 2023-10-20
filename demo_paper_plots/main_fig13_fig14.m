%% Generate Figure 13 & 14 (Hand-Eye Real Data Experiments)
% This code runs for 30-60 minutes.
% Reduce nExp and/or subsetSize for a quicker result.
clc;
clear all;
close all;


%% Load Data
load T_handcam_UR3e210506.mat   % Hand-eye calibration experimental dataset

A = invertData(T_hand);
B = T_cam;


%% Data Normalization
one_rad = (1/pi*180);
scale = 3 * one_rad;

A(1:3,4,:) = A(1:3,4,:) / scale;
B(1:3,4,:) = B(1:3,4,:) / scale;

if size(A,3) ~= size(B,3)
    error('Error loading data.')
end

n = size(A,3);


%% Solving Distance-minimization Using All Data to Acquire Best Solutions
% Parameter Setting
alpha = 2.0;                    % translation weight
param = defaultParam;           % get default solver parameters for distance min.
% param.globalOptMethod = 2;      % activate stochastic global optimization


% Solve with distance-minimization algorithm - conditation 1: AX=YB
[X_geo0{1},Y_geo0{1}] = solveAXYB_sgo(A,B,alpha,param);

% Solve with distance-minimization algorithm - conditation 2: B^-1 Y^-1 = X^-1 A^-1
invA = invertData(A);
invB = invertData(B);
[invY, invX] = solveAXYB_sgo(invB,invA,alpha,param);
X_geo0{2} = invSE3(invX);
Y_geo0{2} = invSE3(invY);

% Solve with distance-minimization algorithm - conditation 3: B X^-1 = Y^-1 A
[invX,invY] = solveAXYB_sgo(B,A,alpha,param);
X_geo0{3} = invSE3(invX);
Y_geo0{3} = invSE3(invY);

% Solve with distance-minimization algorithm - conditation 4: A^-1 Y = X B^-1
[Y_geo0{4},X_geo0{4}] = solveAXYB_sgo(invA,invB,alpha,param);


%% Calculating Error w.r.t. Best Solutions
subsetSize = 10:1:88;
nExp = 500;


X_geo{1} = cell(length(subsetSize), nExp);
Y_geo{1} = cell(length(subsetSize), nExp);

X_geo{2} = cell(length(subsetSize), nExp);
Y_geo{2} = cell(length(subsetSize), nExp);

X_geo{3} = cell(length(subsetSize), nExp);
Y_geo{3} = cell(length(subsetSize), nExp);

X_geo{4} = cell(length(subsetSize), nExp);
Y_geo{4} = cell(length(subsetSize), nExp);

tic;
for j = 1:nExp
    for i = 1:length(subsetSize)
        % random subset of data
        idx = randperm(n);
        idx = idx(1:subsetSize(i));
        A_cur = A(:,:,idx);
        B_cur = B(:,:,idx);
        invA_cur = invA(:,:,idx);
        invB_cur = invB(:,:,idx);
        
        % Solve with distance-minimization algorithm - conf. 1: AX=YB
        [X_geo{1}{i,j},Y_geo{1}{i,j}] = solveAXYB_sgo(A_cur,B_cur,alpha,param);

        % Solve with distance-minimization algorithm - conf. 2: B^-1 Y^-1 = X^-1 A^-1
        [invY, invX] = solveAXYB_sgo(invB_cur,invA_cur,alpha,param);
        X_geo{2}{i,j} = invSE3(invX);
        Y_geo{2}{i,j} = invSE3(invY);

        % Solve with distance-minimization algorithm - conf. 3: B X^-1 = Y^-1 A
        [invX,invY] = solveAXYB_sgo(B_cur,A_cur,alpha,param);
        X_geo{3}{i,j} = invSE3(invX);
        Y_geo{3}{i,j} = invSE3(invY);

        % Solve with distance-minimization algorithm - conf. 4: A^-1 Y = X B^-1
        [Y_geo{4}{i,j},X_geo{4}{i,j}] = solveAXYB_sgo(invA_cur,invB_cur,alpha,param);
    end
    
    t1 = toc;
    disp(['========== ', num2str(j), '/', num2str(nExp), ', time=', num2str(t1), ', time_left=', num2str(t1/j*(nExp-j))]);
end

%% Store Result
for i = 1:length(subsetSize)
    for j = 1:nExp
        % results
        for conf = 1:4
            dX_geo = invSE3(X_geo0{conf}) * X_geo{conf}{i,j};
            dY_geo = invSE3(Y_geo0{conf}) * Y_geo{conf}{i,j};
            distX_geo_SO3{conf}(i,j) = norm(LogSO3(dX_geo(1:3,1:3)));
            distY_geo_SO3{conf}(i,j) = norm(LogSO3(dY_geo(1:3,1:3)));
            distX_geo_trans{conf}(i,j) = norm(dX_geo(1:3,4));
            distY_geo_trans{conf}(i,j) = norm(dY_geo(1:3,4));
        end
        
    end
    
    disp(['============= Storing Results: ', num2str(i), ' out of ', num2str(length(subsetSize))]);
end

%% Plot
clr = lines(4);

figure
subplot(2,1,1)
hold on
for conf = 1:4
    plot(subsetSize, 180/pi*mean(distX_geo_SO3{conf},2), 'Color', clr(conf,:))
end
title('Error in rotation of X')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
ax = gca;
ax.YLim(1) = 0;
xlabel('Number of measurements')
ylabel('Error in rotation of X (^o)')
grid on
% 1-mean(distX_prob_SO3')./mean(distX_geo_SO3')

subplot(2,1,2)
hold on
for conf = 1:4
    plot(subsetSize, scale*mean(distX_geo_trans{conf},2), 'Color', clr(conf,:))
end
title('Error in translation of X')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
ax = gca;
ax.YLim(1) = 0;
xlabel('Number of measurements')
ylabel('Error in translation of X (mm)')
grid on

set(gca,'LooseInset',get(gca,'TightInset'))


figure
subplot(2,1,1)
hold on
for conf = 1:4
    plot(subsetSize, 180/pi*mean(distY_geo_SO3{conf},2), 'Color', clr(conf,:))
end
title('Error in rotation of Y')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
ax = gca;
ax.YLim(1) = 0;
xlabel('Number of measurements')
ylabel('Error in rotation of Y (^o)')
grid on


subplot(2,1,2)
hold on
for conf = 1:4
    plot(subsetSize, scale*mean(distY_geo_trans{conf},2), 'Color', clr(conf,:))
end
title('Error in translation of Y')
legend('Case 1', 'Case 2', 'Case 3', 'Case 4')
ax = gca;
ax.YLim(1) = 0;
xlabel('Number of measurements')
ylabel('Error in translation of Y (mm)')
grid on

set(gca,'LooseInset',get(gca,'TightInset'))

















