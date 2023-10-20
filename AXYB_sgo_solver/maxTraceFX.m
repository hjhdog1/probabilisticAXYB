function X = maxTraceFX( F )
    [U S V] = svd(F);
    detF = det(F);
    
    X = zeros(3);
    if(detF > 0)
        X = V*eye(3)*U';
    else
        X = V*diag([1 1 -1])*U';
    end

end

