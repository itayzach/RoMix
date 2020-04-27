%% Restart run
clear; 
close all; 
clc;
rng('default'); % For reproducibility
%% Parameters
sParams = GetParameters();

%% Clean previous outputs
% ClearPrevOutputs(sParams.sSim.outputFolder);

%% Verify stuff
if sParams.sSim.b_plotEigenFigs    
    [ mPhi_K, vLambda_K ] = CalcAnalyticEigenfunctions(sParams);
    [ mPhi_A, vLambda_A ] = CalcNumericEigenvectors(sParams);
%     PlotEigenfunctionsEigenvectors(sParams, mPhi_K, mPhi_L);
%     PlotSpectrum(sParams, vLambda_AnaL, vLambda_NumL);
    PlotEigenfunctionsEigenvectors(sParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A);
    PlotSpectrum(sParams, vLambda_K, vLambda_A);
end

if sParams.sSim.b_verifyKernelEigenfuncs
    VerifyKernelEigenfuncs(sParams);
else
    warning('Not verifying RKHS');
end

if sParams.sSim.b_verifyEigOrth
    VerifyOrthogonality(sParams);
else
    warning('Not verifying orthogonality');
end

if sParams.sSim.b_verifyMercersTheorem
    VerifyMercerTheorem(sParams);
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
if sParams.sSim.b_extrapolateEnable
    if sParams.dim == 1
        Extrapolate1D(sParams);
    elseif sParams.dim == 2
        Extrapolate2D(sParams);
    else
        error('Cannot extrapolate for more than 2D')
    end
end

