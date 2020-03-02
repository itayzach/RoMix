function [sParams, sSimParams] = GetParameters()

% dimensions
sParams.dim = 1;

% kernel
% sParams.a = 0.1;
% sParams.l = 1/sqrt(2)/10;
sParams.a = 1;
sParams.l = 1/sqrt(2*3); %1/sqrt(2);
sParams.b = 1/(2*sParams.l^2);

% p(x)
sParams.sigma = 1/sqrt(4*sParams.a);
sParams.mu = 0;

% type one consts
sParams.c = sqrt(sParams.a^2 + 2*sParams.a*sParams.b);
sParams.A = sParams.a + sParams.b + sParams.c;
sParams.B = sParams.b/sParams.A;


sParams.omega = 1/sqrt(2*3); % same as sParams.l
sParams.beta = 2*sParams.sigma^2/sParams.omega^2;



% extrapolation
sParams.gamma = 0; % regularization

% num of eigenfunctions
sParams.PlotSpectM = 30;
sParams.RkhsM = 30;
sParams.OrthM = 30;
sParams.MercerM = 100;
sParams.ExtrplM = 10;

% num of sampled points to extrapolate from
sParams.R = 30; %2*M; %floor(N/step);

% simulation
sSimParams.outputFolder = 'figs';
sSimParams.nEigenFuncsToPlot = 3;

sSimParams.b_plotEigenFigs        = true;
sSimParams.b_verifyRKHS           = true;
sSimParams.b_verifyEigOrth        = true;
sSimParams.b_verifyMercersTheorem = true;
sSimParams.b_extrapolateEnable    = true;

sSimParams.b_randomStepSize       = true;

% dataset
sSimParams.twomoons_dataset = false;
sSimParams.twomoons_scale = true;

% AWGN
sSimParams.noiseVar1 = 0;%0.1;
sSimParams.noiseVar2 = 0;%0.1;

assert(sParams.ExtrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
end