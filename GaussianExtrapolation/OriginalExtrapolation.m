%% Restart run
clear; 
close all; 
clc;
rng('default'); % For reproducibility
%% Parameters
sSimParams = GetParameters();

%% Clean previous outputs
% ClearPrevOutputs(sSimParams.outputFolder);

%% Verify stuff
if sSimParams.b_plotEigenFigs    
    [ mPhi_K, vLambda_K ]     = CalcAnalyticEigenfunctions(sSimParams, sSimParams.x_rand);
    [ mPhi_A, vLambda_A ]     = CalcNumericEigenvectors(sSimParams);
    [ mPhi_Nys, vLambda_Nys ] = CalcNystromEigenvectors(sSimParams);
%     PlotEigenfunctionsEigenvectors(sSimParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A);
%     PlotSpectrum(sSimParams, vLambda_K, vLambda_A);
    PlotHistogram(sSimParams);
    PlotEigenDiffs(sSimParams, mPhi_K, mPhi_A, mPhi_Nys)
    PlotEigenfuncvecScatter(sSimParams, mPhi_K, vLambda_K, 'Analytic')
    PlotEigenfuncvecScatter(sSimParams, mPhi_A, vLambda_A, 'Numeric')
    PlotEigenfuncvecScatter(sSimParams, mPhi_Nys, vLambda_Nys, 'Nystrom')
end

if sSimParams.b_verifyKernelEigenfuncs
    VerifyKernelEigenfuncs(sSimParams);
else
    warning('Not verifying RKHS');
end

if sSimParams.b_verifyEigOrth
    VerifyOrthogonality(sSimParams);
else
    warning('Not verifying orthogonality');
end

if sSimParams.b_verifyMercersTheorem
    VerifyMercerTheorem(sSimParams);
else
    warning('Not verifying Mercer''s theorem');
end
%% Extrapolate functions
if sSimParams.b_extrapolateEnable
    if sSimParams.dim == 1
        Extrapolate1D(sSimParams);
    elseif sSimParams.dim == 2
        Extrapolate2D(sSimParams);
    else
        error('Cannot extrapolate for more than 2D')
    end
end

