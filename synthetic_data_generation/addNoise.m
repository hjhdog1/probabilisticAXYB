function R_new = addNoise(R,noiseVar,noisePosition, noiseType)
% R : SO(3)
% noiseVar : noise variance (radian)
% noisePosition : 0-left, 1-right


    if noiseType == 'G'

        if strcmp(noisePosition, 'left')
            R_new = LargeSO3(randn(3,1)*noiseVar)*R;
        elseif strcmp(noisePosition, 'right')
            R_new = R*LargeSO3(randn(3,1)*noiseVar);
        else
            error('Check noise position.');
        end
    elseif noiseType == 'U'
        
        w = randn(3,1);
        w = w/norm(w);
        
        if strcmp(noisePosition, 'left')
            R_new = LargeSO3(w*noiseVar*rand)*R;
        elseif strcmp(noisePosition, 'right')
            R_new = R*LargeSO3(w*noiseVar*rand);
        else
            error('Check noise position.');
        end
    else
        error('Check noise type.');
    end

end

