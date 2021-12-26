function [A,B,M] = generateABData_SE3(X,Y, N, noiseLevel_SO3, position_scaler, noiseLevel_position, outlierRatio, noiseType, noiseAPosition, noiseBPosition)
% N : number of data
% noiseLevel_SO3 : standard deviation for SO3 (radian)
% noiseLevel_position : standard deviation for position
% outlierRation : Ratio of outliers in Data
% M : number of outliers
% noiseAPosition, noiseBPosition: on which side to multiply noise (either 'left' or 'right')

if nargin < 8
    noiseType = 'G';    % Gaussian noise
end

if nargin < 9
    noiseAPosition = 'right';
    noiseBPosition = 'right';
end

M = round(N*outlierRatio);
A = zeros(4,4,N);
B = zeros(4,4,N);
A(4,4,:) = 1;
B(4,4,:) = 1;


isSO3Problem = 0;
if norm(X(1:3,4)) == 0 & norm(Y(1:3,4)) == 0
    isSO3Problem = 0;
end

invX = invSE3(X);
invY = invSE3(Y);

%% Generate noisy data
% 
A0 = randSE3();
A0(1:3,4) = A0(1:3,4)/norm(A0(1:3,4)) * position_scaler * 2.0;

% B0 = randSE3();
% B0(1:3,4) = B0(1:3,4)/norm(B0(1:3,4)) * position_scaler * 2.0;

for i = 1:N-M
    % generate exact A,B
    A(:,:,i) = addNoiseSE3(A0, 0.5, position_scaler, 'right', 'G');
    B(:,:,i) = invY*A(:,:,i)*X;

%     B(:,:,i) = addNoiseSE3(B0, 0.5, position_scaler, 'right', 'G');
%     A(:,:,i) = Y*B(:,:,i)*invX;
    
%     % generate exact A,B
%     A(1:3,1:3,i) = randSO3();
%     if ~isSO3Problem
%         A(1:3,4,i) = position_scaler*randn(3,1);
%     end
%     B(:,:,i) = invY*A(:,:,i)*X;

%     % generate exact A,B
%     B(1:3,1:3,i) = randSO3();
%     if ~isSO3Problem
%         B(1:3,4,i) = position_scaler*randn(3,1);
%     end
%     A(:,:,i) = Y*B(:,:,i)*invX;
    
    % add noise    
    if isSO3Problem
        A(1:3,1:3,i) = addNoise(A(1:3,1:3,i) ,noiseLevel_SO3, noiseAPosition, noiseType);
        B(1:3,1:3,i) = addNoise(B(1:3,1:3,i) ,noiseLevel_SO3, noiseBPosition, noiseType);
    else
        A(:,:,i) = addNoiseSE3(A(:,:,i) ,noiseLevel_SO3, noiseLevel_position, noiseAPosition, noiseType);
        B(:,:,i) = addNoiseSE3(B(:,:,i) ,noiseLevel_SO3, noiseLevel_position, noiseBPosition, noiseType);
    end
%     if ~isSO3Problem
%         invYAX = invY*A(:,:,i)*X;
%         B(1:3,4,i) = invYAX(1:3,4) + noiseLevel_position*randn(3,1);
%         A(1:3,4,i) = A(1:3,4,i) + noiseLevel_position*randn(3,1);
% %         B(1:3,4,i) = B(1:3,4,i) + noiseLevel_position*randn(3,1);
%     end
end


%% Generate outliers
for i = N-M+1:N
    A(:,:,i) = [randSO3() position_scaler*randn(3,1); zeros(1,3) 1];
    B(:,:,i) = [randSO3() position_scaler*randn(3,1); zeros(1,3) 1];
end


end