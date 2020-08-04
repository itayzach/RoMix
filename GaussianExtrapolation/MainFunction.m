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
    [ mPhi_K, vLambda_K ]     = CalcAnalyticEigenfunctions(sParams, sParams.x_rand);
    [ mPhi_A, vLambda_A ]     = CalcNumericEigenvectors(sParams);
    [ mPhi_Nys, vLambda_Nys ] = CalcNystromEigenvectors(sParams);
%     PlotEigenfunctionsEigenvectors(sParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A);
%     PlotSpectrum(sParams, vLambda_K, vLambda_A);
    PlotHistogram(sParams);
    PlotEigenDiffs(sParams, mPhi_K, mPhi_A, mPhi_Nys)
    PlotEigenfuncvecScatter(sParams, mPhi_K, vLambda_K, 'Analytic')
    PlotEigenfuncvecScatter(sParams, mPhi_A, vLambda_A, 'Numeric')
    PlotEigenfuncvecScatter(sParams, mPhi_Nys, vLambda_Nys, 'Nystrom')
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

