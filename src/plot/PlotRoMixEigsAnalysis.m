function PlotRoMixEigsAnalysis(sPlotParams, sPreset, sDataset, sDistParams, sKernelParams, C, Phi, lambdaPhi, PhiInt)
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotNumEigsPerComp
    PlotNumEigsPerComp(sKernelParams);
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc && sPreset.dim <= 3
    plotInd = 0:4;
    [cData{1:numel(plotInd)}] = deal(sDataset.sData.x);
    [cMarkers{1:numel(plotInd)}] = deal('.');
    cSigStr = RepLegend('\\phi', plotInd);
    [cNumCircles{1:numel(plotInd)*2}] = deal((1:sPreset.nLabeled).');
    pltTitle = [];
    PlotGraphSignals(sPlotParams, pltTitle, ...
        ['GaussianKernel_AnalyticEigs_',num2str(plotInd(1)), '_to_', num2str(plotInd(end)) ], cData, ...
        mat2cell(Phi(:,plotInd+1),sPreset.n,ones(1,numel(plotInd))), ...
        cSigStr, cNumCircles, cMarkers, [], [], [min(min(Phi(:,plotInd+1))), max(max(Phi(:,plotInd+1)))]);
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
    vPr = []; %vPr = sDistParams.GMModel.pdf(sDataset.sData.x);
    PlotInnerProductMatrix([], Phi/sqrt(sPreset.n), vPr, '$\frac{1}{\sqrt{n}}{\bf \Phi}^T \frac{1}{\sqrt{n}} {\bf \Phi}$', 'Phi');
    vPr = []; %vPr = sDistParams.GMModel.pdf(sDataset.sData.xt);
    PlotInnerProductMatrix([], PhiInt/sqrt(sPreset.N), vPr, '$\frac{1}{\sqrt{N}}{\bf \Phi}^T \frac{1}{\sqrt{N}} {\bf \Phi}$', 'PhiInt');
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    [~, ~, ~, ~, mAlpha] = InterpGraphSignalRepThm(sPlotParams, sPreset, sDataset);
    CRep = diag(lambdaPhi)*Phi.'*mAlpha;
    PlotCoeffsMatrix(C, '${\bf C}$', CRep, '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
    LnRef = [];
    Cint = RoMix(PhiInt, sPreset.gamma1, sPreset.gamma2, lambdaPhi, LnRef, sDataset.sData.yt);
    PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf int}}$');
    PlotCoefficients(sPlotParams, C, lambdaPhi);
end
end