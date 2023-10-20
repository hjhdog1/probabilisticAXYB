function [A,B,M] = generateABData_SO3(X,Y, N, noiseLevel, outlierRatio)
% N : number of data
% noiseLevel : standard deviation (radian)
% outlierRation : Ratio of outliers in Data
% M : number of outliers

M = round(N*outlierRatio);
A = zeros(3,3,N);
B = zeros(3,3,N);

noiseLevel_radian = noiseLevel;

%% Generate noisy data
for i = 1:N-M
    A(:,:,i) = randSO3();
    B(:,:,i) = Y'*A(:,:,i)*X;
    A(:,:,i) = addNoise(A(:,:,i) ,noiseLevel_radian, 'right', 'G');
    B(:,:,i) = addNoise(B(:,:,i) ,noiseLevel_radian, 'left', 'G');
end


%% Generate outliers
for i = N-M+1:N
    A(:,:,i) = randSO3();
    B(:,:,i) = randSO3();
end


end

