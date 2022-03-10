function [A, dist] = CalcKernel(sKernelParams, x, y)

if ~exist('y', 'var')
    y = x;
end

if strcmp(sKernelParams.kernelType, 'sinc')
    % sinc(x_i - x_j)
    dist = pdist2(x, y);
    A = sin(sKernelParams.a*dist)./(pi*dist);
    A(dist == 0) = sKernelParams.a;
elseif strcmp(sKernelParams.kernelType, 'gaussian')
    dist = pdist2(x, y);
    A = exp(-dist.^2/(2*sKernelParams.omega^2));
else
    error('Unknown kernelType')
end
end