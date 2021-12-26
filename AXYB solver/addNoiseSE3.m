function T_new = addNoiseSE3(T, noiseStd_SO3, noiseStd_position, noisePosition, noiseType)
% T : SE(3)
% noiseStd_SO3 : noise std of rotation matrix (radian)
% noiseStd_position :  : noise std of position (user's length-unit)
% noisePosition : 'left' or 'right'

    if strcmp(noisePosition, 'none')
        T_new = T;
        return;
    end

    T_noise = eye(4);
    if noiseType == 'G'
%         T_noise(1:3, 1:3) = LargeSO3(randn(3,1)*noiseStd_SO3);
%         T_noise(1:3, 1:3) = gaussianRandSO3(eye(3) / (noiseStd_SO3^2));
        T_noise(1:3, 1:3) = gaussianRandSO32(eye(3) / (noiseStd_SO3^2));
        T_noise(1:3, 4) = randn(3,1)*noiseStd_position;
    elseif noiseType == 'U'
        w = randn(3,1);
        w = w/norm(w);        
        T_noise(1:3, 1:3) = LargeSO3(w*noiseStd_SO3);
        
        p = randn(3,1);
        p = p/norm(p);
        T_noise(1:3, 4) = p*noiseStd_position;
    else
        error('Check noise type.');
    end
    
    if strcmp(noisePosition, 'left')
        T_new = T_noise*T;
    elseif strcmp(noisePosition, 'right')
        T_new = T*T_noise;
    else
        error('Check noise position.');
    end

end

