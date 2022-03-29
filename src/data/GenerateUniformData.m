function x = GenerateUniformData(dim, n, xMin, xMax)
x = zeros(n, dim);
for d = 1:dim
    x(:,d) = (xMax(d) - xMin(d))*rand(n, 1) + xMin(d);
end