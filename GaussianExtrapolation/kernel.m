function vKernel = kernel(sKernelParams, x, y) 

if strcmp(sKernelParams.kernelType, 'exp')
    if sKernelParams.constsType == 1
        vKernel = exp(-vecnorm(x-y, 2, 2).^2./(2*sKernelParams.ell^2));
    elseif sKernelParams.constsType == 2
        vKernel = exp(-vecnorm(x-y, 2, 2).^2./(2*sKernelParams.omega^2));
    elseif sKernelParams.constsType == 3
        vKernel = exp(-sKernelParams.eps^2.*vecnorm(x-y, 2, 2).^2);
    else
        error('Unknown constsType');
    end

    % warning('better change kernel to matrix form')

    % n1 = size(x, 1);
    % n2 = size(y, 1);
    % vKernel = exp(-(repmat(sum(x.*x,2)',n2,1) + repmat(sum(y.*y,2),1,n1) - 2*y*x')/(2*l^2));

    % vKernel = prod(vKernel, 2);
elseif strcmp(sKernelParams.kernelType, 'sinc')
%     vKernel = sParams.a^2 - x.*y -0.5*abs(x - y);
    vKernel = sin(sKernelParams.a*(x - y))./(pi*(x - y));
    vKernel((x - y) == 0) = sKernelParams.a/pi;
else
    error('unknown kernelType');
end
end
