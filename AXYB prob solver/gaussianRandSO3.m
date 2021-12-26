function R = gaussianRandSO3(invCov)

while(true)
    R = randSO3;
    w = LogSO3(R);    
    if rand < exp(-0.5*w'*invCov*w)
        return
    end
end

end

