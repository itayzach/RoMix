function [sParams, sSimParams] = GetParameters()

% dimensions
sParams.dim = 2;

% kernel
% sParams.a = 0.1;
% sParams.l = 1/sqrt(2)/10;
sParams.a = 2;
sParams.l = 1/sqrt(2);

% p(x)
sParams.sigma = 1/sqrt(4*sParams.a);
sParams.b = 1/(2*sParams.l^2);

%
sParams.c = sqrt(sParams.a^2 + 2*sParams.a*sParams.b);
sParams.A = sParams.a + sParams.b + sParams.c;
sParams.B = sParams.b/sParams.A;
sParams.lambda_0 = sqrt(2*sParams.a/sParams.A) * sParams.B^0;
sParams.lambda_1 = sqrt(2*sParams.a/sParams.A) * sParams.B^1;


% extrapolation
sParams.gamma = 0; % regularization

% num of eigenfunctions
sParams.M = 150;
sParams.extrplM = 8;

% num of sampled points to extrapolate from
sParams.R = 30; %2*M; %floor(N/step);

% simulation
sSimParams.outputFolder = 'figs';
sSimParams.nEigenFuncsToPlot = 4;

sSimParams.b_plotEigenFigs        = true;
sSimParams.b_verifyRKHS           = false;
sSimParams.b_verifyEigOrth        = false;
sSimParams.b_verifyMercersTheorem = false;
sSimParams.b_extrapolateEnable    = false;

sSimParams.b_randomStepSize       = true;

% dataset
sSimParams.twomoons_dataset = false;
sSimParams.twomoons_scale = true;

% AWGN
sSimParams.noiseVar1 = 0;%0.1;
sSimParams.noiseVar2 = 0;%0.1;

assert(sParams.extrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
end