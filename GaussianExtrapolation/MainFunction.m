%% Restart run
clear; close all; clc;
rng(0); % set seed

%% Parameters
[sParams, sSimParams] = GetParameters();

%% Clean previous outputs
ClearPrevOutputs(sSimParams.outputFolder);

%% Verify stuff
if sSimParams.b_plotEigenFigs
    PlotFirstEigenfunctions(sParams, sSimParams)
end

if sSimParams.b_verifyRKHS
    VerifyRHKS(sParams)
end

if sSimParams.b_verifyEigOrth
    VerifyOrthogonality(sParams)
end

if sSimParams.b_verifyMercersTheorem
    VerifyMercerTheorem(sParams)
end
%% Read mnist data
% mnist = load('data/mnist.mat');
% mXTest = single(mnist.testX);
% vTestLabels = mnist.testY.';
% [nTestPoints, nPixels] = size(mXTest);
% N = nTestPoints;

%% Graph signals
% nDigits = 10;
% mS = zeros(N, nDigits);
% for k = 0:9
%     mS(:, k+1) = (vTestLabels == k);
% end

%% Extrapolate functions
Extrapolate(sParams, sSimParams)


