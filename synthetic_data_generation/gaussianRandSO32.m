function R = gaussianRandSO32(invCov)

A = invCov^(-0.5);

while(true)
    
    w = A * randn(3,1);
    norm_w = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3));
    
    if norm_w > pi
        continue;
    end
    
    if norm_w > 1e-5
        vf = -(2*cos(norm_w) - 2)/(norm_w*norm_w);
    else
        vf = 1.0;
    end
    
    if rand < vf
        R = LargeSO3(w);
        return;
    end
end

end

