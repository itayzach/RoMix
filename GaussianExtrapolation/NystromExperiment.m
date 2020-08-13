%% Restart run
clear; 
close all; 
clc;
rng('default'); % For reproducibility
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureVisible','on')
%% Parameters
vActualDataDist = [ "Gaussian", "Two_moons" ];
vNysRatio = [0.05, 0.4, 0.8];
sSimParams = GetSimParams();
vDataDim = [1 2];
%% Run
T = 10;
nTrain = 5000;
nTest = 0;
for dataDim = vDataDim
    for actualDataDist = vActualDataDist
        if strcmp(actualDataDist, "Two_moons") && dataDim ~= 2
            continue
        end
        for nysRatio = vNysRatio
            tPhiAnalytic = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
            tPhiNumeric  = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
            tPhiNystrom  = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
            for t = 1:T
                sDataset = GenerateDataset(actualDataDist, dataDim, nTrain, nTest);
                sDistParams = EstimateDistributionParameters(sDataset);
                sKernelParams = GetKernelParams(sDataset, sDistParams);
                [sKernelParams.vLambdaAnaytic, sKernelParams.vMultindexToSingleIndexMap] = CalcAnalyticEigenvalues(sSimParams, sKernelParams, sDataset.dim);

                [ tPhiAnalytic(t,:,:), vLambdaAnalytic ]  = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, sDataset.sData.x);
                [ tPhiNumeric(t,:,:), vLambdaNumeric ]    = CalcNumericEigenvectors(sSimParams, sKernelParams, sDataset.sData.x);
                [ tPhiNystrom(t,:,:), vLambdaNystrom ]    = CalcNystromEigenvectors(sSimParams, sKernelParams, sDataset.sData.x, nysRatio);

                [ tPhiNumeric(t,:,:), tPhiNystrom(t,:,:) ] = FlipSign(sSimParams, squeeze(tPhiAnalytic(t,:,:)), squeeze(tPhiNumeric(t,:,:)), squeeze(tPhiNystrom(t,:,:)));

                if t == 1 && sSimParams.b_plotEigenfunctions
                    % Plot histogram and eigenvectors of the first iteration
                    PlotHistogram(sSimParams, sDataset);
                    PlotSpectrum(sSimParams, sDataset, nysRatio, vLambdaAnalytic, vLambdaNumeric, vLambdaNystrom);
                    if dataDim == 1
                        PlotInnerProductMatrix(sSimParams, sDistParams, sDataset, nysRatio, squeeze(tPhiAnalytic(1,:,:)), 'Analytic');
                        PlotInnerProductMatrix(sSimParams, sDistParams, sDataset, nysRatio, squeeze(tPhiNumeric(1,:,:)), 'Numeric');
                        PlotInnerProductMatrix(sSimParams, sDistParams, sDataset, nysRatio, squeeze(tPhiNystrom(1,:,:)), 'Nystrom');
                    end
                    if dataDim == 1
                        firstEigenIdxToPlot = 0;
                        lastEigIdxToPlot = 4;
                        PlotEigenfunctionsEigenvectors(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Analytic')
                        PlotEigenfunctionsEigenvectors(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Nystrom')

                        firstEigenIdxToPlot = 0;
                        lastEigIdxToPlot = sSimParams.PlotEigenFuncsM-1;
                        PlotEigenSamePlot(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), squeeze(tPhiNystrom(1,:,:)))
                        PlotEigenDiffs(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), squeeze(tPhiNystrom(1,:,:)))
                    elseif dataDim == 2
                        firstEigenIdxToPlot = 0;
                        lastEigIdxToPlot = sSimParams.PlotEigenFuncsM-1;
                        PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), vLambdaAnalytic, 'Analytic')
                        PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNumeric(1,:,:)), vLambdaNumeric, 'Numeric')
                        PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(1,:,:)), vLambdaNystrom, 'Nystrom')
                        PlotEigenDiffs2D(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Analytic')
                        PlotEigenDiffs2D(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Analytic')
                    end
                end
            end
            nEig = sSimParams.CalcEigenFuncsM;
            vRMSEAnaVsNum = CalcRMSE(tPhiAnalytic(:,:,1:nEig), tPhiNumeric(:,:,1:nEig), 'Analytic');
            vRMSENysVsNum = CalcRMSE(tPhiNystrom(:,:,1:nEig), tPhiNumeric(:,:,1:nEig), 'Nystrom');        
            PlotRMSE(sSimParams, sDataset, nysRatio, vRMSEAnaVsNum, vRMSENysVsNum);
        end
    end
end