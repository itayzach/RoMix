function vKernel = kernel(x,y,l) 

if ~isvector(x)
    error('x has to be a vector');
end

vKernel = exp(-(x-y).^2./(2*l^2));

% n1 = size(x, 1);
% n2 = size(y, 1);
% vKernel = exp(-(repmat(sum(x.*x,2)',n2,1) + repmat(sum(y.*y,2),1,n1) - 2*y*x')/(2*l^2));

% vKernel = prod(vKernel, 2);
end
