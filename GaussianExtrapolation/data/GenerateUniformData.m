function x = GenerateUniformData(dim, n)
xMin = -1;
xMax = 1;
x = (xMax - xMin)*rand(n, dim) + xMin;
end