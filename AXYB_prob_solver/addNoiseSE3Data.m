function [A, B] = addNoiseSE3Data(A, B, noiseLevel_SO3, noiseLevel_position, noiseAPosition, noiseBPosition, noiseType)

n = size(A,3);
for i = 1:n
    A(:,:,i) = addNoiseSE3(A(:,:,i) ,noiseLevel_SO3, noiseLevel_position, noiseAPosition, noiseType);
    B(:,:,i) = addNoiseSE3(B(:,:,i) ,noiseLevel_SO3, noiseLevel_position, noiseBPosition, noiseType);
end

end

