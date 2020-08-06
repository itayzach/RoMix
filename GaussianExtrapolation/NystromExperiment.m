%% Restart run
clear; 
close all; 
clc;
rng('default'); % For reproducibility
%% Parameters
vActualDataDist = [ "Gaussian", "Two_moons" ];
vNysRatio = [0.8, 0.05, 0.4];
sSimParams = GetSimParams();

%% Run
T = 50;
nTrain = 5000;
nTest = 0;
for actualDataDist = vActualDataDist
    for nysRatio = vNysRatio
        tPhiAnalytic = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
        tPhiNumeric  = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
        tPhiNystrom  = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
        for t = 1:T
            sDataset = GenerateDataset(actualDataDist, nTrain, nTest);
            sDistParams = EstimateDistributionParameters(sDataset);
            sKernelParams = GetKernelParamsAndCalcEigenvalues(sSimParams, sDataset, sDistParams);
            
            [ tPhiAnalytic(t,:,:), vLambdaAnalytic ]  = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, sDataset.mData.x);
            [ tPhiNumeric(t,:,:), vLambdaNumeric ]    = CalcNumericEigenvectors(sSimParams, sKernelParams, sDataset.mData.x);
            [ tPhiNystrom(t,:,:), vLambdaNystrom ]    = CalcNystromEigenvectors(sSimParams, sKernelParams, sDataset.mData.x, nysRatio);
            
            [ tPhiNumeric, tPhiNystrom ] = FlipSign(sSimParams, tPhiAnalytic, tPhiNumeric, tPhiNystrom);
            
            if t == 1 && sSimParams.b_plotEigenfunctions
                % Plot histogram and eigenvectors of the first iteration
                PlotHistogram(sSimParams, sDataset);
%                 PlotInnerProductMatrix(sSimParams, squeeze(tPhiAnalytic(1,:,:)), 'Analytic');
                PlotInnerProductMatrix(sSimParams, sDataset, nysRatio, squeeze(tPhiNumeric(1,:,:)), 'Numeric');
                PlotInnerProductMatrix(sSimParams, sDataset, nysRatio, squeeze(tPhiNystrom(1,:,:)), 'Nystrom');
                firsfirstEigenIdxToPlot = 0;
                lastEigIdxToPlot = 19;
                PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), vLambdaAnalytic, 'Analytic')
                PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNumeric(1,:,:)), vLambdaNumeric, 'Numeric')
                PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(1,:,:)), vLambdaNystrom, 'Nystrom')
                
                if strcmp(actualDataDist, 'Two_moons') && nysRatio == 0.05
                    firsfirstEigenIdxToPlot = 30;
                    lastEigIdxToPlot = 49;
                    PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), vLambdaAnalytic, 'Analytic')
                    PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNumeric(1,:,:)), vLambdaNumeric, 'Numeric')
                    PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firsfirstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(1,:,:)), vLambdaNystrom, 'Nystrom')
                end
            end
        end
        
        if strcmp(actualDataDist, 'Two_moons') && nysRatio == 0.05
            vRMSEAnaVsNum = CalcRMSE(tPhiAnalytic(:,:,1:10), tPhiNumeric(:,:,1:10), 'Analytic');
            vRMSENysVsNum = CalcRMSE(tPhiNystrom(:,:,1:10), tPhiNumeric(:,:,1:10), 'Nystrom');
        else
            vRMSEAnaVsNum = CalcRMSE(tPhiAnalytic, tPhiNumeric, 'Analytic');
            vRMSENysVsNum = CalcRMSE(tPhiNystrom, tPhiNumeric, 'Nystrom');
        end
        PlotEigenDiffs(sSimParams, sDataset, nysRatio, vRMSEAnaVsNum, vRMSENysVsNum);
        
    end
end
