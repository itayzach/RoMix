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
        
        pos=find(sParams.sDataset.y==1);
        neg=find(sParams.sDataset.y==-1);
        sParams.sDataset.x(pos,:) = sParams.sDataset.x(pos,:) + [ 0 0.25 ];
        sParams.sDataset.x(neg,:) = sParams.sDataset.x(neg,:) - [ 0 0.25 ];
        
        pos=find(sParams.sDataset.yt==1);
        neg=find(sParams.sDataset.yt==-1);
        sParams.sDataset.xt(pos,:) = sParams.sDataset.xt(pos,:) + [ 0 0.25 ];
        sParams.sDataset.xt(neg,:) = sParams.sDataset.xt(neg,:) - [ 0 0.25 ];
    end
end
%% p(x)
sParams.dataDist = 'gaussian';

if strcmp(sParams.dataDist, 'gaussian')
    if isfield(sParams, 'sDataset')
        GMModel = fitgmdist(sParams.sDataset.x,1);
        sParams.cov = GMModel.Sigma;
        sParams.mu  = GMModel.mu;
    else
        if sParams.dim == 2
            sParams.cov = [0.25    0; 
                           0   0.25];
            sParams.mu  = [0 0];
        elseif sParams.dim == 1
            sParams.cov = 0.5;
            sParams.mu = 100;
        else
            error('generate random data for more than 2D')
        end
    end
    [sParams.u, sParams.sigma_eigv] = eig(sParams.cov);
%     [sParams.u, sParams.sigma_eigv] = svd(sParams.cov);
    sParams.sigma = diag(sqrt(sParams.sigma_eigv)).';
    if sParams.dim == 2
        if sParams.cov(1,1) > sParams.cov(2,2)
            warning('The variance in the first axis is greater than the variance in the second, but eig returns the eigenvalues in increasing order. So we fliplr')
            sParams.u = fliplr(sParams.u);    
            sParams.sigma = fliplr(sParams.sigma);
        end   
        sParams.u = [-sParams.u(:,1) -sParams.u(:,2)];    
    end
    sParams.mu_1D = sParams.mu*sParams.u;
    isalmostequal(sParams.u*diag(sParams.sigma.^2)*sParams.u.', sParams.cov, 1e-15)

    sParams.xMax = sParams.mu + 1.5*sParams.sigma;
    sParams.xMin = sParams.mu - 1.5*sParams.sigma;
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
% if sParams.sSim.twomoons_dataset
%     x_rand = [sParams.sDataset.x; sParams.sDataset.xt];
% else
    n = 5000;
    if strcmp(sParams.dataDist, 'gaussian')
        x_rand = mvnrnd(sParams.mu, sParams.cov, n);
    elseif strcmp(sParams.dataDist, 'uniform')
        x_rand = (sParams.xMax - sParams.xMin)*rand(n, sParams.dim) + sParams.xMin;
    else
        error('unknown pdf')
    end
% end
sParams.x_rand = x_rand;

%% x-axis
sParams.n_x_axis = 100;
sParams.dx = (sParams.xMax - sParams.xMin)/sParams.n_x_axis;
x = zeros(sParams.n_x_axis, sParams.dim);
for d = 1:sParams.dim
    x(:,d) = (sParams.xMin(d) : sParams.dx(d) : sParams.xMax(d)-sParams.dx(d)).';
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

        sParams.omega = 0.3;1/(6*sqrt(2)); % kernel width
        sParams.beta = 2*sParams.sigma.^2/sParams.omega^2;
        sParams.t = 0.5*sParams.omega^2;

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
sParams.PlotEigenFuncsM = 23;
sParams.PlotSpectM = 30;
sParams.RkhsM = 20;
sParams.OrthM = 30;
sParams.MercerM = 0;

sParams.ExtrplM = 20;

sParams.FirstM = 0;
%% extrapolation
sParams.gamma = 0; % regularization
sParams.R = 5000;    % num of sampled points to extrapolate from

%% simulation
sParams.sSim.outputFolder = 'figs';

sParams.sSim.b_plotEigenFigs          = false;
sParams.sSim.b_verifyKernelEigenfuncs = false;
sParams.sSim.b_verifyEigOrth          = false;
sParams.sSim.b_verifyMercersTheorem   = false;
sParams.sSim.b_extrapolateEnable      = false;

sParams.sSim.b_randomStepSize = true;
sParams.sSim.b_plot_contourf = false;

%% AWGN
sParams.sSim.noiseVar1 = 0; %0.1;
sParams.sSim.noiseVar2 = 0; %0.1;

%% Get correct order of eigenvalues (for 1D indexing from multindexing)
if sParams.dim == 2
    vLambda_K = zeros(max(sParams.PlotSpectM,sParams.ExtrplM),1);
    for i = 0:max(sParams.PlotSpectM,sParams.ExtrplM)-1
        m = OneDim2TwoDimIndex(i);
        vLambda_K(i+1) = lambda(sParams, m);
    end
end
[sParams.vLambda_K, sParams.multindexToSingleIndexMap] = sort(vLambda_K, 'descend');


end