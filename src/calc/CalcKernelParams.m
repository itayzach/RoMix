function [sKernelParams, t] = CalcKernelParams(sDistParams, omega)
ts = tic;
%% kernel and eigenfunctions constants type
if strcmp(sDistParams.estDataDist, 'Gaussian')
    sKernelParams.kernelType = 'gaussian';
    %% second type consts
    sKernelParams.omega = omega;
    sKernelParams.t = 0.5*sKernelParams.omega^2;
    for c = 1:sDistParams.estNumComponents
        sKernelParams.beta{c} = 2*sDistParams.sigma{c}.^2/sKernelParams.omega^2;
    end
    PrintKernelParams(sDistParams,sKernelParams);
elseif strcmp(sDistParams.estDataDist, 'uniform')
    sKernelParams.kernelType = 'sinc';
else
    error('unknown pdf')
end

%% Save distParams
sKernelParams.sDistParams = sDistParams;
t = toc(ts);

end