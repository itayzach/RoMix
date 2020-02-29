%% Restart run
clear; 
close all; 
clc;
rng(0); % set seed

%% Parameters
[sParams, sSimParams] = GetParameters();

%% Clean previous outputs
% ClearPrevOutputs(sSimParams.outputFolder);

%% Verify stuff
if sSimParams.b_plotEigenFigs    

    [ mPhi_A, vLambda_A ] = PlotNumericEigenvectors(sParams, sSimParams);
    [ mPhi_K, vLambda_K ] = PlotAnalyticEigenfunctions(sParams, sSimParams);
    PlotSpectrum(sParams, sSimParams, vLambda_A);

end

if sSimParams.b_verifyRKHS
    VerifyRHKS(sParams);
else
    warning('Not verifying RKHS');
end

if sSimParams.b_verifyEigOrth
    VerifyOrthogonality(sParams);
else
    warning('Not verifying orthogonality');
end

if sSimParams.b_verifyMercersTheorem
    VerifyMercerTheorem(sParams, sSimParams);
else
    warning('Not verifying Mercer''s theorem');
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
if sSimParams.b_extrapolateEnable
    if sParams.dim == 1
        Extrapolate1D(sParams, sSimParams);
    elseif sParams.dim == 2
        Extrapolate2D(sParams, sSimParams);
    else
        error('Cannot extrapolate for more than 2D')
    end
end

