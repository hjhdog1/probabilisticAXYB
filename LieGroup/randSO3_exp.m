function R = randSO3_exp()
% generate random SO(3)

while true
    while true
        w = 2*pi * (rand(3,1)-0.5);
        norm_w = norm(w);
        if(norm_w < pi)
            break;
        end 
    end
    
    if (norm_w < 1e-12)
        R = eye(3);
        return;
    end
    
    p = (2-2*cos(norm_w))/(norm_w^2);
    
    if p > rand
        R = LargeSO3(w);
        return;
    end

end


end

