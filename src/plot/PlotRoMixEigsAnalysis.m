function PlotRoMixEigsAnalysis(sPlotParams, sPreset, sDataset, sDistParams, sKernelParams, C, Phi, lambdaPhi, PhiInt)
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotNumEigsPerComp
    PlotNumEigsPerComp(sKernelParams);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc && sPreset.dim <= 3
    plotInd = [0 4]*(sPreset.dim == 1) + [0 11]*(sPreset.dim > 1);
    figTitle = 'Eigenfunctions of the Gaussian kernel (on $n$ nodes)';
    figName = 'Phi';
    PlotKernelEigenfunctions([], 'Gaussian', sDataset.sData.x, [], plotInd(1), plotInd(end), Phi, C, 'C', [], figTitle, figName, '\phi' );
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc
    figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    PlotSpectrum([], [], lambdaPhi, [], [], '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotMercer
    W = CalcAdjacency(sDataset.sData.x, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
    CheckMercerTheorem(sDistParams, Phi, lambdaPhi, W, sDataset.sData.x);
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
    Cint = RoMix(PhiInt, sPreset.gamma1, sPreset.gamma2, lambdaPhi, LnRef, sDataset.sData.yt);
    PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf int}}$');
end
end