function sKernelParams = GetKernelParams(sDistParams, omega)
dim = sDistParams.dim;
%% kernel and eigenfunctions constants type
if strcmp(sDistParams.estDataDist, 'Gaussian')
    sKernelParams.kernelType = 'gaussian';
    %% second type consts
    fprintf('*            Using beta,omega constants                 *\n');
    fprintf('*********************************************************\n');
    sKernelParams.omega = omega;
    sKernelParams.t = 0.5*sKernelParams.omega^2;
    for c = 1:sDistParams.estNumComponents
        sKernelParams.beta{c} = 2*sDistParams.sigma{c}.^2/sKernelParams.omega^2;

        for d = 1:dim
            fprintf('sigma(%d) (pdf width)     = %8.3f\n', d, sDistParams.sigma{c}(d));
            fprintf('mu(%d)    (pdf mean)      = %8.3f\n', d, sDistParams.mu{c}(d));
        end
            fprintf('omega    (kernel width)  = %8.3f\n', sKernelParams.omega);
        for d = 1:dim
            fprintf('--> beta(%d)              = %8.3f\n', d, sKernelParams.beta{c}(d));
        end
    end
    fprintf('*********************************************************\n');
elseif strcmp(sDistParams.estDataDist, 'uniform')
    sKernelParams.kernelType = 'sinc';
else
    error('unknown pdf')
end

%% Save distParams
sKernelParams.sDistParams = sDistParams;


end