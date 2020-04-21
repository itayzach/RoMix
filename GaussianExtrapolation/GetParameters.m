function [sParams] = GetParameters()

%% dimensions
sParams.dim = 2;
fprintf('*********************************************************\n');
fprintf('*                      %d-D                              *\n', sParams.dim);
fprintf('*********************************************************\n');
%% Dataset
if sParams.dim == 2
    sParams.sSim.twomoons_dataset = true;
    if sParams.sSim.twomoons_dataset
        sParams.sDataset = load('2moons.mat');
        sParams.sSim.twomoons_scale = false;
    end
end
%% p(x)
sParams.dataDist = 'gaussian';

if strcmp(sParams.dataDist, 'gaussian')
    if isfield(sParams, 'sDataset')
        GMModel = fitgmdist([sParams.sDataset.x; sParams.sDataset.xt],1);
        sParams.cov = GMModel.Sigma;
        [sParams.u, sParams.sigma_eigv] = eig(sParams.cov);
        if sParams.cov(1,1) > sParams.cov(2,2)
            warning('The variance in the first axis is greater than the variance in the second, but eig returns the eigenvalues in increasing order. So we fliplr')
            sParams.u = fliplr(sParams.u);    
            sParams.sigma = [sParams.sigma_eigv(2,2) sParams.sigma_eigv(1,1)];
        else
            sParams.sigma = [sParams.sigma_eigv(1,1) sParams.sigma_eigv(2,2)];
        end
        sParams.mu = GMModel.mu;
    else
        sParams.sigma = 0.5*ones(1, sParams.dim);
        sParams.cov   = diag(sParams.sigma);
        sParams.mu    = 0*ones(1, sParams.dim);
        sParams.u     = eye(sParams.dim);
    end
    
    sParams.xMax = sParams.mu + 3*sParams.sigma;
    sParams.xMin = sParams.mu - 3*sParams.sigma;
    if sParams.xMax - sParams.xMin > 10
        sParams.xMax = sParams.mu + 5;
        sParams.xMin = sParams.mu - 5;
        warning('3sigma is too large...');
    end
elseif strcmp(sParams.dataDist, 'uniform')
    sParams.a    = 2.5;
    sParams.u    = eye(sParams.dim);
    sParams.xMax = sParams.a;
    sParams.xMin = -sParams.a;
else
    error('unknown pdf')
end
%% random x-axis
n = 5000;
if strcmp(sParams.dataDist, 'gaussian')
    x_rand = sParams.sigma.*randn(n, sParams.dim) + sParams.mu;
    sParams.xMin = max(sParams.xMin, min(x_rand));
    sParams.xMax = min(sParams.xMax, max(x_rand));
    warning('Note: sParams.xMin, sParams.xMax are changed');
elseif strcmp(sParams.dataDist, 'uniform')
    x_rand = (sParams.xMax - sParams.xMin)*rand(n, sParams.dim) + sParams.xMin;
else
    error('unknown pdf')
end
sParams.x_rand = x_rand;

%% x-axis
sParams.n_x_axis = 1000;
sParams.dx = (sParams.xMax(1) - sParams.xMin(1))/sParams.n_x_axis;

x = zeros(sParams.n_x_axis, sParams.dim);
for d = 1:sParams.dim
    x(:,d) = (sParams.xMin : sParams.dx : sParams.xMax-sParams.dx).';
end
sParams.x = x;

%% kernel and eigenfunctions constants type
if strcmp(sParams.dataDist, 'gaussian')
    sParams.kernelType = 'exp';
    sParams.constsType = 2;
    if sParams.constsType == 1
        %% first type consts
        fprintf('*              Using a,b,ell constants                  *\n');
        fprintf('*********************************************************\n');

        sParams.ell = 1/sqrt(2); % kernel width

        sParams.a = 1./(2*sParams.sigma);
        sParams.b = 1/(2*sParams.ell^2);

        sParams.c = sqrt(sParams.a.^2 + 2*sParams.a.*sParams.b);
        sParams.A = sParams.a + sParams.b + sParams.c;
        sParams.B = sParams.b./sParams.A;

        for d = 1:sParams.dim
            fprintf('a(%d) = %8.3f --> sigma(%d) (pdf width)    = %8.3f\n', d, sParams.a(d), d, sParams.sigma(d));
        end
        fprintf('b     = %8.3f --> ell   (kernel width) = %8.3f\n', sParams.b, sParams.ell); 
        fprintf('*********************************************************\n');
    elseif sParams.constsType == 2
        %% second type consts
        fprintf('*            Using beta,omega constants                 *\n');
        fprintf('*********************************************************\n');

        sParams.omega = 1/(6*sqrt(2)); % kernel width
        sParams.beta = 2*sParams.sigma.^2/sParams.omega^2;

        for d = 1:sParams.dim
            fprintf('sigma(%d) (pdf width)     = %8.3f\n', d, sParams.sigma(d));
            fprintf('mu(%d)    (pdf mean)      = %8.3f\n', d, sParams.mu(d));
        end
            fprintf('omega    (kernel width)  = %8.3f\n', sParams.omega);
        for d = 1:sParams.dim
            fprintf('--> beta(%d)              = %8.3f\n', d, sParams.beta(d));
        end
        fprintf('*********************************************************\n');
    elseif sParams.constsType == 3
        %% third type consts
        fprintf('*            Using alpha,eps constants                  *\n');
        fprintf('*********************************************************\n');

        sParams.eps = 1; % 1/kernel width
        sParams.alpha = 1./(sqrt(2)*sParams.sigma);

        for d = 1:sParams.dim
            fprintf('alpha(d)                  = %8.3f\n', d, sParams.alpha(d));
        end
        fprintf('eps (kernel width)        = %8.3f\n', sParams.eps);
        for d = 1:sParams.dim
            fprintf('--> sigma(%d) (pdf width)  = %8.3f\n', d, sParams.sigma(d));
            fprintf('    mu(%d)    (pdf mean)   = %8.3f\n', d, sParams.mu(d));
        end
        fprintf('*********************************************************\n');    
    else
        error('Unknown constsType')
    end
elseif strcmp(sParams.dataDist, 'uniform')
    sParams.kernelType = 'sinc';
else
    error('unknown pdf')
end


%% num of eigenfunctions
sParams.PlotEigenFuncsM = 4;
sParams.PlotSpectM = 30;
sParams.RkhsM = 20;
sParams.OrthM = 30;
sParams.MercerM = 2500;
sParams.ExtrplM = 2500;
sParams.FirstM = 2000;

%% extrapolation
sParams.gamma = 0; % regularization
sParams.R = 5000;    % num of sampled points to extrapolate from

%% simulation
sParams.sSim.outputFolder = 'figs';

sParams.sSim.b_plotEigenFigs          = false;
sParams.sSim.b_verifyKernelEigenfuncs = true;
sParams.sSim.b_verifyEigOrth          = false;
sParams.sSim.b_verifyMercersTheorem   = false;
sParams.sSim.b_extrapolateEnable      = false;

sParams.sSim.b_randomStepSize       = true;


%% AWGN
sParams.sSim.noiseVar1 = 0; %0.1;
sParams.sSim.noiseVar2 = 0; %0.1;

assert(sParams.ExtrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
end