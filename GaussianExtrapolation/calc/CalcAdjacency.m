function [A, dist] = CalcAdjacency(sKernelParams, x, y)

if ~exist('y', 'var')
    y = x;
end

if strcmp(sKernelParams.kernelType, 'sinc')
    % sinc(x_i - x_j)
    dist = pdist2(x, y);
    A = sin(sKernelParams.a*dist)./(pi*dist);
    A(dist == 0) = sKernelParams.a;
elseif strcmp(sKernelParams.kernelType, 'gaussian')
    % exp(-||x_i - x_j||^2 / 2\ell^2)
    dist = pdist2(x, y);
    if sKernelParams.constsType == 1
        A = exp(-dist.^2/(2*sKernelParams.ell^2));
    elseif sKernelParams.constsType == 2
        A = exp(-dist.^2/(2*sKernelParams.omega^2));
    elseif sKernelParams.constsType == 3
        A = exp(-sKernelParams.eps^2 * dist.^2);
    else
        error('Unknown constsType');
    end
else
    error('Unknown kernelType')
end
end