function PlotEigenvectorsWrapper(sPlotParams, sPreset, sDataset, lambdaPhi, Phi, PhiInt, LnRef, VRefToCompare, VIntToCompare, C)
xInt = sDataset.sData.xt;
% ----------------------------------------------------------------------------------------------
% Plots
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    b_maskDataTermCMatrix = false;
    Cint = RoMix(PhiInt, sPreset.gamma1, sPreset.gamma2, lambdaPhi, LnRef, VRefToCompare, b_maskDataTermCMatrix);
    PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf int}}$');
end
if sPlotParams.b_globalPlotEnable && sPreset.dim <= 3
    if sPreset.dim == 1
        plotInd = 0:4;
        PlotEigenfuncvecScatter(sPlotParams, sPreset.verticesPDF, xInt, [], plotInd(1), plotInd(end), ...
            VRefToCompare, [], [], [], ['Eigenvectors interpolation (N = ', num2str(sPreset.N), ')'], ...
            'VInt', '\tilde{\psi}', VIntToCompare, '\tilde{\psi}^{{\bf RoMix}}');
    else
        plotInd = 1:4;
        [cData{1:numel(plotInd)*2}] = deal(xInt);
        cSigStr = [RepLegend('\\tilde{\\psi}', plotInd), RepLegend('\\tilde{\\psi}^{{\\bf RoMix}}', plotInd)];
        [cNumCircles{1:numel(plotInd)*2}] = deal((1:sPreset.N).');
        [cMarkers{1:numel(plotInd)*2}] = deal('.');
        PlotGraphSignals(sPlotParams, ['Eigenvectors interpolation (N = ', num2str(sPreset.N), ')'], ...
            [sPreset.matrixForEigs, '_Eigs_',num2str(plotInd(1)), '_to_', num2str(plotInd(end)) ], cData, ...
            [mat2cell(VRefToCompare(:,plotInd+1),sPreset.N,ones(1,numel(plotInd))), mat2cell(VIntToCompare(:,plotInd+1),sPreset.N,ones(1,numel(plotInd)))], ...
            cSigStr, cNumCircles, cMarkers, [], [], [min(min(VRefToCompare(:,plotInd+1))), max(max(VRefToCompare(:,plotInd+1)))]);
    end
end
end