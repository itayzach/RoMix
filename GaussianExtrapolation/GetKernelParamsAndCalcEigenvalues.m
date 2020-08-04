function sKernelParams = GetKernelParamsAndCalcEigenvalues(sSimParams, sDataset, sDistParams)
%% kernel and eigenfunctions constants type
if strcmp(sDataset.estDataDist, 'Gaussian')
    sKernelParams.kernelType = 'exp';
    sKernelParams.constsType = 2;
    if sKernelParams.constsType == 1
        %% first type consts
        fprintf('*              Using a,b,ell constants                  *\n');
        fprintf('*********************************************************\n');

        sKernelParams.ell = 1/sqrt(2); % kernel width

        sKernelParams.a = 1./(2*sDistParams.sigma);
        sKernelParams.b = 1/(2*sKernelParams.ell^2);

        sKernelParams.c = sqrt(sKernelParams.a.^2 + 2*sKernelParams.a.*sKernelParams.b);
        sKernelParams.A = sKernelParams.a + sKernelParams.b + sKernelParams.c;
        sKernelParams.B = sKernelParams.b./sKernelParams.A;

        for d = 1:sDataset.dim
            fprintf('a(%d) = %8.3f --> sigma(%d) (pdf width)    = %8.3f\n', d, sKernelParams.a(d), d, sDistParams.sigma(d));
        end
        fprintf('b     = %8.3f --> ell   (kernel width) = %8.3f\n', sKernelParams.b, sKernelParams.ell); 
        fprintf('*********************************************************\n');
    elseif sKernelParams.constsType == 2
        %% second type consts
        fprintf('*            Using beta,omega constants                 *\n');
        fprintf('*********************************************************\n');

        sKernelParams.omega = 0.3; 1/(6*sqrt(2)); % kernel width
        sKernelParams.beta = 2*sDistParams.sigma.^2/sKernelParams.omega^2;
        sKernelParams.t = 0.5*sKernelParams.omega^2;

        for d = 1:sDataset.dim
            fprintf('sigma(%d) (pdf width)     = %8.3f\n', d, sDistParams.sigma(d));
            fprintf('mu(%d)    (pdf mean)      = %8.3f\n', d, sDistParams.mu(d));
        end
            fprintf('omega    (kernel width)  = %8.3f\n', sKernelParams.omega);
        for d = 1:sDataset.dim
            fprintf('--> beta(%d)              = %8.3f\n', d, sKernelParams.beta(d));
        end
        fprintf('*********************************************************\n');
    elseif sKernelParams.constsType == 3
        %% third type consts
        fprintf('*            Using alpha,eps constants                  *\n');
        fprintf('*********************************************************\n');

        sKernelParams.eps = 1; % 1/kernel width
        sKernelParams.alpha = 1./(sqrt(2)*sDistParams.sigma);

        for d = 1:sDataset.dim
            fprintf('alpha(d)                  = %8.3f\n', d, sKernelParams.alpha(d));
        end
        fprintf('eps (kernel width)        = %8.3f\n', sKernelParams.eps);
        for d = 1:sDataset.dim
            fprintf('--> sigma(%d) (pdf width)  = %8.3f\n', d, sDistParams.sigma(d));
            fprintf('    mu(%d)    (pdf mean)   = %8.3f\n', d, sDistParams.mu(d));
        end
        fprintf('*********************************************************\n');    
    else
        error('Unknown constsType')
    end
elseif strcmp(sDataset.estDataDist, 'uniform')
    sKernelParams.kernelType = 'sinc';
else
    error('unknown pdf')
end

%% Save distParams
sKernelParams.sDistParams = sDistParams;
%% Calculate eigenvalues
[sKernelParams.vLambdaAnaytic, sKernelParams.vMultindexToSingleIndexMap] = CalcAnalyticEigenvalues(sSimParams, sKernelParams, sDataset.dim);


end