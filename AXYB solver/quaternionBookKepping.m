function signs = quaternionBookKepping(A)

N = size(A,3);
q = zeros(4,N);
signs = ones(1,N);

q(:,1) = SO3toQuart(A(:,:,1));
for i = 2:N
    q(:,i) = SO3toQuart(A(:,:,i));
    
    positiveDiff = norm(q(:,i) - q(:,i-1));
    negativeDiff = norm(q(:,i) - (-q(:,i-1)));
    
    if positiveDiff < negativeDiff
        signs(i) = signs(i-1);
    else
        signs(i) = -signs(i-1);
    end
end

end

