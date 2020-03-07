function [sParams, sSimParams] = GetParameters()

%% dimensions
sParams.dim = 1;

%% kernel and eigenfunctions constants type
sParams.kernelType = 'exp'; % exp only. sinc isn't yet supported
sParams.constsType = 1;

if sParams.constsType == 1
    %% first type consts
    fprintf('*********************************************************\n');
    fprintf('*              Using a,b,ell constants                  *\n');
    fprintf('*********************************************************\n');

    sParams.a = 1;
    sParams.b = 3;
    
    % kernel width
    sParams.ell = 1/sqrt(2*sParams.b);
    
    % p(x)
    sParams.sigma = 1/(2*sParams.a)*ones(1, sParams.dim);
    sParams.mu = 0*ones(1, sParams.dim);

    sParams.c = sqrt(sParams.a^2 + 2*sParams.a*sParams.b);
    sParams.A = sParams.a + sParams.b + sParams.c;
    sParams.B = sParams.b/sParams.A;
    
    sParams.xMax = sParams.mu + 3*sParams.sigma;
    sParams.xMin = sParams.mu - 3*sParams.sigma;
    fprintf('a = %8.3f --> sigma (pdf width)    = %8.3f\n', sParams.a, sParams.sigma);
    fprintf('b = %8.3f --> ell   (kernel width) = %8.3f\n', sParams.b, sParams.ell); 
    fprintf('*********************************************************\n');
elseif sParams.constsType == 2
    %% second type consts
    fprintf('*********************************************************\n');
    fprintf('*            Using beta,omega constants                 *\n');
    fprintf('*********************************************************\n');
    
    % p(x)
    sParams.sigma = 1*ones(1, sParams.dim);
    sParams.mu    = 0*ones(1, sParams.dim);

    % sParams.sigma = [0.9788    0.4815];
    % sParams.mu    = [0.6858    0.2503];
    
    sParams.xMax = sParams.mu + 3*sParams.sigma;
    sParams.xMin = sParams.mu - 3*sParams.sigma;
    
    sParams.omega = 1/sqrt(2); % kernel width
    sParams.beta = 2*sParams.sigma.^2/sParams.omega^2;
    fprintf('sigma (pdf width)    = %8.3f\n', sParams.sigma);
    fprintf('omega (kernel width) = %8.3f\n', sParams.omega);
    fprintf('--> beta             = %8.3f\n', sParams.beta);
    fprintf('*********************************************************\n');
else
    error('Unknown constsType')
end

%% pdf
sParams.pdf = 'gaussian';

%% num of eigenfunctions
sParams.PlotEigenFuncsM = 4;
sParams.PlotSpectM = 30;
sParams.RkhsM = 20;
sParams.OrthM = 30;
sParams.MercerM = 50;
sParams.ExtrplM = 10;

%% extrapolation
sParams.gamma = 0; % regularization
sParams.R = 30;    % num of sampled points to extrapolate from

%% simulation
sSimParams.outputFolder = 'figs';

sSimParams.b_plotEigenFigs        = true;
sSimParams.b_verifyRKHS           = false;
sSimParams.b_verifyEigOrth        = false;
sSimParams.b_verifyMercersTheorem = false;
sSimParams.b_extrapolateEnable    = true;

sSimParams.b_randomStepSize       = true;

%% dataset
sSimParams.twomoons_dataset = true;
sSimParams.twomoons_scale = true;

%% AWGN
sSimParams.noiseVar1 = 0;%0.1;
sSimParams.noiseVar2 = 0;%0.1;

assert(sParams.ExtrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
end