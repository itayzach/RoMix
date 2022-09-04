function PlotRoMixWrapper(sPlotParams, sPreset, sDataset, sDistParams, sKernelParams, Phi, C, lambdaPhi, mSigCnvrtRec, mSigCnvrtInt, PhiInt)
clusterMethod = 'GMM';
xTrain = sDataset.sData.x;
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotDataVsGmm 
    nGmmPoints = 50*ismember(sPreset.verticesPDF, {'USPS', 'MNIST'}) + sPreset.n*~ismember(sPreset.verticesPDF, {'USPS', 'MNIST'});
    pltTitle = ['Dataset with n = ', num2str(sPreset.n), ' points'];
    if strcmp(clusterMethod, 'GMM')
        plt2Title = ['Generated ' num2str(nGmmPoints), ' points from ' clusterMethod ' with nEstComp = ' num2str(sPreset.gmmNumComponents)];
    elseif strcmp(clusterMethod, 'SC')
        plt2Title = ['Painted ' num2str(sPreset.n), ' according to ' clusterMethod ' with nEstComp = ' num2str(sPreset.gmmNumComponents)];
    end
    windowStyle = 'normal';
    PlotDataset(sPlotParams, sPreset, xTrain, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle);
    if ismember(sPreset.verticesPDF, {'SwissRoll'}) && sPreset.sDatasetParams.b_randn
        b_transform = true;
        PlotDataset(sPlotParams, sPreset, sDataset.sData.S, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, b_transform);
    end
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc && sPreset.dim <= 3
    plotInd = [0 4]*(sPreset.dim == 1) + [0 11]*(sPreset.dim > 1);
    figTitle = 'Eigenfunctions of the Gaussian kernel (on $n$ nodes)';
    figName = 'Phi';
    PlotEigenfuncvecScatter([], 'Gaussian', xTrain, [], plotInd(1), plotInd(end), Phi, C, 'C', [], figTitle, figName, '\phi' );
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc
    figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    PlotSpectrum([], [], lambdaPhi, [], [], '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotClustersAnalysis && sPreset.dim > 1
    PlotGaussianEllipses(sPlotParams, sDistParams);
    PlotCovEigs(sPlotParams, sDistParams);
    if sDistParams.GMModel.NumComponents < 100
        PlotClustersMeans(sPreset, sPlotParams, sDistParams);
    end
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotNumEigsPerComp
    PlotNumEigsPerComp(sKernelParams);
end
if sPlotParams.b_globalPlotEnable
    mSigCnvrtRecRef = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigCnvrtIntRef = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    if sPreset.b_runGraphSignals
        PlotGraphSignalsWrapper(sPlotParams, sPreset, sKernelParams, sDataset, ...
            mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt, C, 'RoMix')
    end
    if sPreset.b_interpEigenvecs
        PlotEigenvectorsWrapper(sPlotParams, sPreset, sDataset, mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt);
    end
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotMercer
    W = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
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
    LnRef = [];
    Cint = RoMix(PhiInt, sPreset.gamma1, sPreset.gamma2, lambdaPhi, LnRef, mSigCnvrtIntRef);
    PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf int}}$');
end
end