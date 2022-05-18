function [mSigCnvrtRec, mSigCnvrtInt, tTrain, tInt, C] = InterpGraphSignalRoMix(sPlotParams, sPreset, sDataset)
xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end

% 1. Estimate distribution parameters (GMM)
[sDistParams, tTrainVec(1:2)] = EstimateDistributionParameters(xTrain, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);
if isfield(sDataset, 'vae'), sDistParams.vae = sDataset.vae; end

% 2. Calculate Phi(xTrain) and lambdaPhi
[sKernelParams, tTrainVec(3)] = CalcKernelParams(sDistParams, sPreset.omegaTilde);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex, tTrainVec(4)] = CalcAnalyticEigenvalues(sPreset.MTilde, sKernelParams);
[ Phi, lambdaPhi ] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xTrain);

% 3. Find C
Ln = [];
[C, tTrainVec(4)] = RoMix(Phi, sPreset.gamma1, sPreset.gamma2, lambdaPhi, Ln, mSig, sPreset.b_maskDataFitTerm);

% 4. Project (reconstruct)
mSigRecPhi = Phi*C;
mSigCnvrtRec = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecPhi);

% 5. Calculate Phi(xInt) and interpolate
[PhiInt, ~, tIntVec(1)] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xInt);
ts = tic;
mSigInt = PhiInt*C;
mSigCnvrtInt = ConvertSignalByDataset(sPreset.verticesPDF, mSigInt);
tIntVec(2) = toc(ts);

tTrain = sum(tTrainVec);
tInt = sum(tIntVec);
fprintf('RoMix: train took %.2f sec, interp took %.2f sec\n', tTrain, tInt);

% ----------------------------------------------------------------------------------------------
% Plots
% ----------------------------------------------------------------------------------------------
clusterMethod = 'GMM';
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotDataVsGmm 
    nGmmPoints = 50*ismember(sPreset.verticesPDF, {'USPS', 'MNIST'}) + sPreset.n*~ismember(sPreset.verticesPDF, {'USPS', 'MNIST'});
    pltTitle = ['Dataset with n = ', num2str(sPreset.n), ' points'];
    if strcmp(clusterMethod, 'GMM')
        plt2Title = ['Generated ' num2str(nGmmPoints), ' points from ' clusterMethod ' with nEstComp = ' num2str(sPreset.gmmNumComponents)];
    elseif strcmp(clusterMethod, 'SC')
        plt2Title = ['Painted ' num2str(sPreset.n), ' according to ' clusterMethod ' with nEstComp = ' num2str(sPreset.gmmNumComponents)];
    end
    windowStyle = 'normal';
    PlotDataset(sPlotParams, xTrain, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc && sPreset.dim <= 3
    plotInd = [0 4]*(sPreset.dim == 1) + [0 11]*(sPreset.dim > 1);
    figTitle = 'Eigenfunctions of the Gaussian kernel (on $n$ nodes)';
    figName = 'Phi';
    PlotEigenfuncvecScatter([], 'Gaussian', xTrain, [], plotInd(1), plotInd(end), Phi, [], [], [], figTitle, figName, '\phi' );
    figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    PlotSpectrum([], [], lambdaPhi, [], [], '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotClustersAnalysis && sPreset.dim > 1
    PlotGaussianEllipses(sPlotParams, sDistParams);
    PlotCovEigs(sPlotParams, sDistParams);
    if sDistParams.GMModel.NumComponents < 100
        PlotClustersMeans(sPreset, sDistParams);
    end
end
if sPlotParams.b_globalPlotEnable
    mSigCnvrtRecRef = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigCnvrtIntRef = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    if sPreset.b_runGraphSignals
        PlotGraphSignalsWrapper(sPlotParams, sPreset, sKernelParams, sDataset, ...
            mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt, C, 'RoMix')
    end
    if sPreset.b_interpEigenvecs
        LnRef = [];
        PlotEigenvectorsWrapper(sPlotParams, sPreset, sDataset, lambdaPhi, Phi, PhiInt, LnRef, mSigCnvrtIntRef, mSigCnvrtInt, C)
    end
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotMercer
    [W, tW, dist, D, Ln, tLn] = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
    CheckMercerTheorem(sDistParams, Phi, lambdaPhi, W, xTrain);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotInnerProductMatrices
    PlotInnerProductMatrix([], Phi/sqrt(sPreset.n), [], '$\frac{1}{\sqrt{n}}{\bf \Phi}^T \frac{1}{\sqrt{n}} {\bf \Phi}$', 'Phi');
    PlotInnerProductMatrix([], PhiInt/sqrt(sPreset.N), [], '$\frac{1}{\sqrt{N}}{\bf \Phi}^T \frac{1}{\sqrt{N}} {\bf \Phi}$', 'PhiInt');
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    [~, ~, ~, ~, mAlpha] = InterpGraphSignalRepThm(sPlotParams, sPreset, sDataset);
    CRep = diag(lambdaPhi)*Phi.'*mAlpha;
    PlotCoeffsMatrix(C, '${\bf C}$', CRep, '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
end
end

