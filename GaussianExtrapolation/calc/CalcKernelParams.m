function sKernelParams = CalcKernelParams(sDistParams, omega)
dim = sDistParams.dim;
%% kernel and eigenfunctions constants type
if strcmp(sDistParams.estDataDist, 'Gaussian')
    sKernelParams.kernelType = 'gaussian';
    %% second type consts
    fprintf('*********************************************************\n');
    fprintf('*            Using beta,omega constants                 *\n');
    fprintf('*********************************************************\n');
    sKernelParams.omega = omega;
    sKernelParams.t = 0.5*sKernelParams.omega^2;
    fprintf('omega    (kernel width)            = %8.3f\n', sKernelParams.omega);
    fprintf('---------------------------------------------------------\n');
    for c = 1:sDistParams.estNumComponents
        sKernelParams.beta{c} = 2*sDistParams.sigma{c}.^2/sKernelParams.omega^2;
        fprintf('mu(%d)    (pdf mean)                = [',c);
        for d = 1:dim
            fprintf('%8.3f, ', sDistParams.mu{c}(d)); 
        end
        fprintf(']\n');
        fprintf('sigma(%d) (pdf width)               = [',c);
        for d = 1:dim
            fprintf('%8.3f, ', sDistParams.sigma{c}(d)); 
        end
        fprintf(']\n');
        fprintf('--> beta(%d) = 2*sigma(%d)^2/omega^2 = [', c, c);
        for d = 1:dim
            fprintf('%8.3f, ', sKernelParams.beta{c}(d));
        end
        fprintf(']\n');
        fprintf('---------------------------------------------------------\n');
    end
elseif strcmp(sDistParams.estDataDist, 'uniform')
    sKernelParams.kernelType = 'sinc';
else
    error('unknown pdf')
end

%% Save distParams
sKernelParams.sDistParams = sDistParams;


end