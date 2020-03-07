function vKernel = kernel(sParams, x, y) 

if ~isvector(x)
    error('x has to be a vector');
end
if sParams.constsType == 1
    vKernel = exp(-(x-y).^2./(2*sParams.l^2));
elseif sParams.constsType == 2
    vKernel = exp(-(x-y).^2./(2*sParams.omega^2));
else
    error('Unknown constsType');
end
% n1 = size(x, 1);
% n2 = size(y, 1);
% vKernel = exp(-(repmat(sum(x.*x,2)',n2,1) + repmat(sum(y.*y,2),1,n1) - 2*y*x')/(2*l^2));

% vKernel = prod(vKernel, 2);
end
