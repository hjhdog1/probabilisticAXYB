function inverted_cov_tensor = invert_cov(cov_tensor)
% This functions inverts all matrics in cov_tensor

inverted_cov_tensor = zeros(size(cov_tensor));

n = size(cov_tensor,3);
for i = 1:n
    inverted_cov_tensor(:,:,i) = inv(cov_tensor(:,:,i));
end

end