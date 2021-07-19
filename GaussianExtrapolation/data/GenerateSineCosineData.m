function [x, theta] = GenerateSineCosineData(dim, n)

theta = rand(n,1);
assert(mod(dim,2) == 0)
x = zeros(n,dim);
for k = 1:dim/2
    x(:,2*k-1) = cos(2*pi*k*theta);
    x(:,2*k) = sin(2*pi*k*theta);
end
end