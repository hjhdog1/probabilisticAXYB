function pairs = allPossiblePair(set)

n = length(set);
pairs = zeros(round(n*(n-1)/2), 2);

count = 1;
for i = 1:n
    for j = i+1:n
        pairs(count, :) = [set(i), set(j)];
        count = count + 1;
    end
end

end

