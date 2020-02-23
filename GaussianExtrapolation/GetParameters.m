function [sParams, sSimParams] = GetParameters()

% dimensions
sParams.dim = 1;

% kernel
sParams.a = 1;
sParams.b = 3;

% p(x)
sParams.sigma = 1/sqrt(4*sParams.a);
sParams.l = 1/sqrt(2*sParams.b);

% extrapolation
sParams.gamma = 0; % regularization

% num of eigenfunctions
sParams.M = 10;

% num of sampled points to extrapolate from
sParams.R = 100; %2*M; %floor(N/step);

% simulation
sSimParams.outputFolder = 'figs';
sSimParams.nEigenFuncsToPlot = 6;

sSimParams.b_plotEigenFigs        = true;
sSimParams.b_verifyRKHS           = false;
sSimParams.b_verifyEigOrth        = false;
sSimParams.b_verifyMercersTheorem = false;
sSimParams.b_extrapolateEnable    = false;

sSimParams.b_randomStepSize       = true;

% AWGN
sSimParams.noiseVar1 = 0;%0.1;
sSimParams.noiseVar2 = 0;%0.1;
end