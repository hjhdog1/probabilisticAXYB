function T = randSE3(position_scaler)
% noiseLevel_position : standard deviation for random positon
    
    if ~exist('position_scaler')
        position_scaler = 1;
    end

    p = position_scaler*randn(3,1);
    T = [randSO3() p; zeros(1,3) 1];

end

