function Rs = randMultipleSO3(N)
% N : number of SO3 to be generated
    Rs = zeros(3,3,N);
    
    for i = 1:N
        Rs(:,:,i) = randSO3;
    end


end

