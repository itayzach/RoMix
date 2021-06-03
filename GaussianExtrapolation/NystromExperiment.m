%% Restart run
clear; close all; clc;
rng('default'); % For reproducibility
set(0,'DefaultFigureWindowStyle','normal')
set(0,'DefaultFigureVisible','on')
%% Parameters
interpMethod = 'NewPoints';
vActualDataDist = [ "Gaussian", "TwoMoons" ];

vNysRatio = [0.8 0.4 0.05];
actualNumComponents = 2;
estNumComponents = actualNumComponents;
sSimParams = GetSimParams();
vDataDim = [1 2];

%% Run
b_normalize = true;
T = 10;
nTrain = 10000;
nTest = 0;
nTotal = nTrain + nTest;
R = length(vNysRatio);
for dataDim = vDataDim
    for actualDataDist = vActualDataDist
        if strcmp(actualDataDist, "TwoMoons") && dataDim ~= 2
            continue
        end
        tPhiAnalytic = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
        tPhiNumeric  = zeros(T, nTrain, sSimParams.CalcEigenFuncsM);
        tPhiNystrom  = zeros(R, T, nTrain, sSimParams.CalcEigenFuncsM);
        for t = 1:T
            if actualNumComponents == 2 && strcmp(actualDataDist, "Gaussian")
                sDatasetParams.mu = [2; 7];       % Means
                sDatasetParams.sigma = [0.4 0.5]; % Covariances
            else
                sDatasetParams = [];
            end
            sDataset = GenerateDataset(actualDataDist, dataDim, actualNumComponents, nTrain, nTest, interpMethod, sDatasetParams);
            if strcmp(actualDataDist, "TwoMoons")
                GMMRegVal = 0.1;
            else
                GMMRegVal = 0;
            end
            sDistParams = EstimateDistributionParameters(sDataset, estNumComponents, GMMRegVal);
            sKernelParams = GetKernelParams(sDataset, sDistParams);
            [sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
                = CalcAnalyticEigenvalues(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.dim, estNumComponents);
            [ tPhiAnalytic(t,:,:), vLambdaAnalytic ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.x, b_normalize);
            [ tPhiNumeric(t,:,:), vLambdaNumeric ]   = CalcNumericEigenvectors(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.x);
            tPhiNumeric(t,:,:) = FlipSign(squeeze(tPhiAnalytic(t,:,:)), squeeze(tPhiNumeric(t,:,:)));
            mLambdaNystrom = zeros(R, sSimParams.CalcEigenFuncsM);
            for r = 1:R
                nysRatio = vNysRatio(r);
                [ tPhiNystrom(r,t,:,:), mLambdaNystrom(r,:) ]   = CalcNystromEigenvectors(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.x, nysRatio); 
                tPhiNystrom(r,t,:,:) = FlipSign(squeeze(tPhiAnalytic(t,:,:)), squeeze(tPhiNystrom(r,t,:,:)));
            end

            if t == 1 && sSimParams.b_plotEigenfunctions
                % Plot histogram and eigenvectors of the first iteration
                PlotHistogram(sSimParams, [ sDataset.sData.x; sDataset.sData.xt], actualDataDist, 'Histogram');
                PlotSpectrum(sSimParams, sDataset, vNysRatio, vLambdaAnalytic, vLambdaNumeric, ...
                    mLambdaNystrom, '\lambda^{{\bf ana}}', '\lambda^{{\bf num}}', '\lambda^{{\bf nys}}');

                graphName = 'ipmatrix';
                PlotInnerProductMatrix(sSimParams, sDistParams.dim, sDistParams.vPr, graphName, [], squeeze(tPhiAnalytic(1,:,:)), 'Analytic')
                PlotInnerProductMatrix(sSimParams, sDistParams.dim, sDistParams.vPr, graphName, [], squeeze(tPhiNumeric(1,:,:)), 'Numeric')
                for r = 1:R
                    nysRatio = vNysRatio(r);
                    PlotInnerProductMatrix(sSimParams, sDistParams.dim, sDistParams.vPr, graphName, nysRatio, squeeze(tPhiNystrom(r,1,:,:)), 'Nystrom')
                end
                if dataDim == 1
                    firstEigenIdxToPlot = 0;
                    lastEigIdxToPlot = 4;
                    PlotEigenfunctionsEigenvectors(sSimParams, sDataset, [], firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Analytic')
                    for r = 1:R
                        nysRatio = vNysRatio(r);
                        PlotEigenfunctionsEigenvectors(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(r,1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Nystrom')
                    end
                    firstEigenIdxToPlot = 0;
                    lastEigIdxToPlot = sSimParams.PlotEigenFuncsM-1;
                    b_plotErrVsNodeInd = true;
                    PlotEigenDiffs(sSimParams, sDataset.dataDist, sDataset.sData.x, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Phi_vs_V', '\phi', 'v', b_plotErrVsNodeInd)
                    for r = 1:R
                        nysRatio = vNysRatio(r);
                        PlotEigenSamePlot(sSimParams, sDataset, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), squeeze(tPhiNystrom(r,1,:,:)))
                        
                        PlotEigenDiffs(sSimParams, sDataset.dataDist, sDataset.sData.x, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(r,1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Nys_vs_V', '\hat{v}', 'v', b_plotErrVsNodeInd)
                    end
                elseif dataDim == 2
                    firstEigenIdxToPlot = 0;
                    lastEigIdxToPlot = sSimParams.PlotEigenFuncsM-1;
                    b_plotErrVsNodeInd = true;
                    PlotEigenfuncvecScatter(sSimParams, sDataset.actualDataDist, sDataset.sData.x, [], firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), vLambdaAnalytic, [], [], 'Eigenfunctions', 'Phi', '\phi')
                    PlotEigenfuncvecScatter(sSimParams, sDataset.actualDataDist, sDataset.sData.x, [], firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNumeric(1,:,:)), vLambdaNumeric, [], [], 'Eigenvectors', 'V', 'v')
                    PlotEigenDiffs(sSimParams, sDataset.dataDist, sDataset.sData.x, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiAnalytic(1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Phi_vs_V', '\phi', 'v', b_plotErrVsNodeInd)
                    for r = 1:R
                        nysRatio = vNysRatio(r);
                        PlotEigenfuncvecScatter(sSimParams, sDataset.actualDataDist, sDataset.sData.x, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(r,1,:,:)), mLambdaNystrom(r,:), [], [], 'Eigenvectors', 'Nys', '\hat{v}')
                        PlotEigenDiffs(sSimParams, sDataset.dataDist, sDataset.sData.x, nysRatio, firstEigenIdxToPlot, lastEigIdxToPlot, squeeze(tPhiNystrom(r,1,:,:)), squeeze(tPhiNumeric(1,:,:)), 'Nys_vs_V', '\hat{v}', 'v', b_plotErrVsNodeInd)
                    end
                end
            end
        end
        nEig = sSimParams.CalcEigenFuncsM;
        vRMSEAnaVsNum = CalcRMSE(tPhiAnalytic(:,:,1:nEig), tPhiNumeric(:,:,1:nEig), 'Analytic');
        mRMSENysVsNum = zeros(R, sSimParams.CalcEigenFuncsM);
        for r = 1:R
            nysRatio = vNysRatio(r);
            mRMSENysVsNum(r,:) = CalcRMSE(reshape(tPhiNystrom(r,:,:,1:nEig),[T, nTotal, nEig]), tPhiNumeric(:,:,1:nEig), 'Nystrom');        
        end
        PlotRMSE(sSimParams, sDataset, vNysRatio, vRMSEAnaVsNum, mRMSENysVsNum);
    end
end