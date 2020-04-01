function vKernel = kernel(sParams, x, y) 

if sParams.constsType == 1
    vKernel = exp(-vecnorm(x-y, 2, 2).^2./(2*sParams.ell^2));
elseif sParams.constsType == 2
    vKernel = exp(-vecnorm(x-y, 2, 2).^2./(2*sParams.omega^2));
elseif sParams.constsType == 3
    vKernel = exp(-sParams.eps^2.*vecnorm(x-y, 2, 2).^2);
else
    error('Unknown constsType');
end

warning('better change kernel to matrix form')

% n1 = size(x, 1);
% n2 = size(y, 1);
% vKernel = exp(-(repmat(sum(x.*x,2)',n2,1) + repmat(sum(y.*y,2),1,n1) - 2*y*x')/(2*l^2));

% vKernel = prod(vKernel, 2);
end
