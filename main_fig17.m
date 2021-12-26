%% Generate Figure 13 & 14 (Dual-camera Real Data Experiments)
% This code runs for 10-20 hours.
% Reduce nExp for a quicker result.
clc;
clear all;
close all;


%% Load Data
load T_b1_c1.mat
A = T;

load T_b2_c2.mat
B = T;

%% Normalize Data
one_rad = (1/pi*180);
scale = 3 * one_rad;

A(1:3,4,:) = A(1:3,4,:) / scale;
B(1:3,4,:) = B(1:3,4,:) / scale;

if size(A,3) ~= size(B,3)
    error('Error loading data.')
end

n = size(A,3);


%% Solve with distance-minimization algorithm
% Parameter Setting
alpha = 2.0;                    % translation weight
param = defaultParam;           % get default solver parameters for distance min.
% param.globalOptMethod = 2;      % activate stochastic global optimization

% Solve AX = YB with geometric stochastic global optimization
[X_geo0,Y_geo0] = solveAXYB_SE3(A,B,alpha,param);

X_true = X_geo0;
Y_true = Y_geo0;
    
%% Solve with probabilistic algorithm
invSig_wN = repmat(eye(3), [1,1,n]);
invSig_pN = repmat(eye(3), [1,1,n]);
invSig_wM = repmat(eye(3), [1,1,n]);
invSig_pM = repmat(eye(3), [1,1,n]);

noiseConf = 2;

step_R = 0.01;
step_p = 0.01;

[X_prob0, Y_prob0, C, ~, N, M] = solveAXYB_prob(A, B, X_geo0, Y_geo0, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p);


%% Calculating Error w.r.t. Geometric Solution
subsetSize = 10:5:100;
nExp = 500;


X_geo = cell(length(subsetSize), nExp);
Y_geo = cell(length(subsetSize), nExp);

X_prob = cell(length(subsetSize), nExp);
Y_prob = cell(length(subsetSize), nExp);

tic;
for j = 1:nExp
    for i = 1:length(subsetSize)
        % random subset of data
        idx = randperm(n);
        idx = idx(1:subsetSize(i));
        A_cur = A(:,:,idx);
        B_cur = B(:,:,idx);
        
        % solve AX=YB
        [X_geo{i,j}, Y_geo{i,j}] = solveAXYB_SE3(A_cur, B_cur ,alpha ,param);
        [X_prob{i,j}, Y_prob{i,j}, C] = solveAXYB_prob(A_cur, B_cur, X_geo{i,j}, Y_geo{i,j}, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p);
    end
    
    t1 = toc;
    disp(['========== ', num2str(j), '/', num2str(nExp), ', time=', num2str(t1), ', time_left=', num2str(t1/j*(nExp-j))]);
end


%% Result
% Process Polaris Data
folder = './data/dualCam/NDI_polaris/data210504_1/';

p = cell(2,2);
q = cell(2,2);

for i = 1:2
    for j = 1:2
        tbl_p = readtable([folder, 'p', num2str(i-1), num2str(j-1), '.csv']);
        tbl_q = readtable([folder, 'q', num2str(i-1), num2str(j-1), '.csv']);
        
        p{i,j} = mean(tbl_p{:,10:12});
        q{i,j} = mean(tbl_q{:,10:12});
        
    end
end

p = cell2mat(p(:))';
q = cell2mat(q(:))';

dpq = reshape(repmat(p,4,1),3,[]) - repmat(q,1,4);
d = sqrt(sum((dpq).*(dpq)));

w = 12.4571;
for i = 1:length(subsetSize)
    for j = 1:nExp
        % results of geo
        p_geo = w * [0, 0, 5, 5; 0, 4, 0, 4; 0, 0, 0, 0];
        X_cur = X_geo{i,j};
        X_cur(1:3,4) = X_cur(1:3,4) * scale;
        q_geo = X_cur * [p_geo; ones(1,4)];
        q_geo = q_geo(1:3,:);
        
        dpq_geo = reshape(repmat(p_geo,4,1),3,[]) - repmat(q_geo,1,4);
        d_geo = sqrt(sum((dpq_geo).*(dpq_geo)));
        distX_geo(i,j) = mean(abs(d_geo - d));
        
        % results of prob
        p_prob = w * [0, 0, 5, 5; 0, 4, 0, 4; 0, 0, 0, 0];
        X_cur = X_prob{i,j};
        X_cur(1:3,4) = X_cur(1:3,4) * scale;
        q_prob = X_cur * [p_prob; ones(1,4)];
        q_prob = q_prob(1:3,:);
        
        dpq_prob = reshape(repmat(p_prob,4,1),3,[]) - repmat(q_prob,1,4);
        d_prob = sqrt(sum((dpq_prob).*(dpq_prob)));
        distX_prob(i,j) = mean(abs(d_prob - d));

        
    end
    
    disp(['============= Storing Results: ', num2str(i), ' out of ', num2str(length(subsetSize))]);
end


%% Plot
clr = lines(2);

figure('Position', [100,100,500,200])
plot(subsetSize, mean(distX_geo,2), 'Color', clr(1,:))
hold on
plot(subsetSize, mean(distX_prob,2), 'Color', clr(2,:))
grid on
% title('Average error in distance between corner points of two checkerboards')
xlabel('Number of measurements')
ylabel('Average distance error (mm)')
legend('Distance minimization', 'Proposed method')
ax = gca;
ax.YLim(1) = 0;

set(gca,'LooseInset',get(gca,'TightInset'))

