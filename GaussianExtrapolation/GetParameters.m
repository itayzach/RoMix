function [sParams, sSimParams] = GetParameters()

%% dimensions
sParams.dim = 2;
fprintf('*********************************************************\n');
fprintf('*                      %d-D                              *\n', sParams.dim);
fprintf('*********************************************************\n');
%% kernel and eigenfunctions constants type
sParams.kernelType = 'exp'; % exp only. sinc isn't yet supported
sParams.constsType = 1;

if sParams.constsType == 1
    %% first type consts
    fprintf('*              Using a,b,ell constants                  *\n');
    fprintf('*********************************************************\n');

    % kernel width
    sParams.ell = 1/sqrt(2)/10;
%     sParams.ell = 1/sqrt(2*3);
    sParams.b = 1/(2*sParams.ell^2);
    
    % p(x)
    sParams.sigma = 500*ones(1, sParams.dim);
%     sParams.sigma = 0.5*ones(1, sParams.dim);
    sParams.mu = 0*ones(1, sParams.dim);
    sParams.a = 1./(2*sParams.sigma);

    sParams.c = sqrt(sParams.a.^2 + 2*sParams.a.*sParams.b);
    sParams.A = sParams.a + sParams.b + sParams.c;
    sParams.B = sParams.b./sParams.A;
    
    fprintf('a = %8.3f --> sigma (pdf width)    = %8.3f\n', sParams.a, sParams.sigma);
    fprintf('b = %8.3f --> ell   (kernel width) = %8.3f\n', sParams.b, sParams.ell); 
    fprintf('*********************************************************\n');
elseif sParams.constsType == 2
    %% second type consts
    fprintf('*            Using beta,omega constants                 *\n');
    fprintf('*********************************************************\n');
    
    % p(x)
    sParams.sigma = 0.5*ones(1, sParams.dim);
    sParams.mu    = 0*ones(1, sParams.dim);

    % sParams.sigma = [0.9788    0.4815];
    % sParams.mu    = [0.6858    0.2503];
        
    sParams.omega = 1/sqrt(2); % kernel width
    sParams.beta = 2*sParams.sigma.^2/sParams.omega^2;
    fprintf('sigma (pdf width)    = %8.3f\n', sParams.sigma);
    fprintf('omega (kernel width) = %8.3f\n', sParams.omega);
    fprintf('--> beta             = %8.3f\n', sParams.beta);
    fprintf('*********************************************************\n');
elseif sParams.constsType == 3
    %% third type consts
    fprintf('*            Using alpha,eps constants                  *\n');
    fprintf('*********************************************************\n');
    
    sParams.eps = 1; % 1/kernel width
    sParams.alpha = 2/sqrt(2)*ones(1, sParams.dim);
    
    % p(x)
    sParams.sigma = 1./(sqrt(2)*sParams.alpha);
    sParams.mu    = 0*ones(1, sParams.dim);

    % sParams.sigma = [0.9788    0.4815];
    % sParams.mu    = [0.6858    0.2503];
        
    fprintf('alpha     = %8.3f\n', sParams.alpha);
    fprintf('epsilon   = %8.3f\n', sParams.eps);
    fprintf('--> sigma (pdf width)  = %8.3f\n', sParams.sigma);
    fprintf('*********************************************************\n');    
else
    error('Unknown constsType')
end

%% pdf
sParams.dataDist = 'gaussian';
sParams.xMax = sParams.mu + 3*sParams.sigma;
sParams.xMin = sParams.mu - 3*sParams.sigma;
if sParams.xMax - sParams.xMin > 10
    sParams.xMax = sParams.mu + 5;
    sParams.xMin = sParams.mu - 5;
    warning('3sigma is too large...');
end
%% x-axis
sParams.n_x_axis = 1000;
sParams.dx = (sParams.xMax(1) - sParams.xMin(1))/sParams.n_x_axis;

x = zeros(sParams.n_x_axis, sParams.dim);
for d = 1:sParams.dim
    x(:,d) = (sParams.xMin : sParams.dx : sParams.xMax-sParams.dx).';
end
sParams.x = x;
%% random x-axis
n = 5000;
x_rand = zeros(n, sParams.dim);
if strcmp(sParams.dataDist, 'gaussian')
    x_rand = sParams.sigma.*randn(n, sParams.dim) + sParams.mu;
    sParams.xMin = max(sParams.xMin, min(x_rand));
    sParams.xMax = min(sParams.xMax, max(x_rand));
    warning('Note: sParams.xMin, sParams.xMax are changed');
elseif strcmp(sParams.dataDist, 'uniform')
    sParams.xMin = -0.5;
    sParams.xMax = 0.5;
    x_rand(:, d) = (sParams.xMax - sParams.xMin)*rand(n, 1) + sParams.xMin;
    warning('Note: sParams.xMin, sParams.xMax are changed');
else
    error('unknown pdf')
end
sParams.x_rand = x_rand;



%% num of eigenfunctions
sParams.PlotEigenFuncsM = 9;
sParams.PlotSpectM = 30;
sParams.RkhsM = 20;
sParams.OrthM = 30;
sParams.MercerM = 50;
sParams.ExtrplM = 10;

%% extrapolation
sParams.gamma = 0; % regularization
sParams.R = 50;    % num of sampled points to extrapolate from

%% simulation
sSimParams.outputFolder = 'figs';

sSimParams.b_plotEigenFigs        = true;
sSimParams.b_verifyRKHS           = false;
sSimParams.b_verifyEigOrth        = true;
sSimParams.b_verifyMercersTheorem = false;
sSimParams.b_extrapolateEnable    = true;

sSimParams.b_randomStepSize       = true;

%% dataset
sSimParams.twomoons_dataset = false;
sSimParams.twomoons_scale = true;

%% AWGN
sSimParams.noiseVar1 = 0; %0.1;
sSimParams.noiseVar2 = 0; %0.1;

assert(sParams.ExtrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
end