%% Restart run
clear; 
close all; 
clc;
rng('default'); % For reproducibility
%% Parameters
vNysRatio = [0.05, 0.4, 0.8];
sSimParams = GetSimParams();

%% Run
T = 50;
nTrain = 5000;
nTest = 0;
for simIdx = 1:2
    if simIdx == 1
        actualDataDist = 'Gaussian';
    else
        actualDataDist = 'Two_moons';
    end
    for nysRatio = vNysRatio
        tPhiAnalytic = zeros(T, nTrain, sSimParams.PlotEigenFuncsM);
        tPhiNumeric  = zeros(T, nTrain, sSimParams.PlotEigenFuncsM);
        tPhiNystrom  = zeros(T, nTrain, sSimParams.PlotEigenFuncsM);
        for t = 1:T
            sDataset = GenerateDataset(actualDataDist, nTrain, nTest);
            sDistParams = EstimateDistributionParameters(sDataset);
            sKernelParams = GetKernelParamsAndCalcEigenvalues(sSimParams, sDataset, sDistParams);
            
            [ tPhiAnalytic(t,:,:), vLambdaAnalytic ]  = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, sDataset.mData.x);
            [ tPhiNumeric(t,:,:), vLambdaNumeric ]    = CalcNumericEigenvectors(sSimParams, sKernelParams, sDataset.mData.x);
            [ tPhiNystrom(t,:,:), vLambdaNystrom ]    = CalcNystromEigenvectors(sSimParams, sKernelParams, sDataset.mData.x, nysRatio);
            
            [ tPhiNumeric, tPhiNystrom ] = FlipSign(sSimParams, tPhiAnalytic, tPhiNumeric, tPhiNystrom);
            
            if t == 1
                % Plot histogram and eigenvectors of the first iteration
                PlotHistogram(sSimParams, sDataset)
                firsfirstEigenIdxToPlot = 0;
                lastEigIdxToPlot = 19;
                PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), vLambdaAnalytic, 'Analytic')
                PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNumeric(1,:,:)), vLambdaNumeric, 'Numeric')
                PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(1,:,:)), vLambdaNystrom, 'Nystrom')
            end
        end
        
        vRMSEAnaVsNum = CalcRMSE(sSimParams, T, tPhiAnalytic, tPhiNumeric);
        vRMSENysVsNum = CalcRMSE(sSimParams, T, tPhiNystrom, tPhiNumeric);
        PlotEigenDiffs(sSimParams, sDataset, nysRatio, vRMSEAnaVsNum, vRMSENysVsNum);
        
    end
end