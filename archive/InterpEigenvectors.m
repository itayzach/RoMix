function [VIntToCompare, VRecPhi, VToCompare, VNysToCompare, VRepToCompare, VRefToCompare] = ...
    InterpEigenvectors(sPlotParams, sPreset, sDataset, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, VRef, W, WRef, D, DRef, Ln, LnRef)
xInt = sDataset.sData.xt;
interpRatio = sPreset.N/sPreset.n;
% ----------------------------------------------------------------------------------------------
% Interpolate with our method
% ----------------------------------------------------------------------------------------------
if sPreset.b_forceCtoIdentity
    C = zeros(sPreset.MTilde, sPreset.M);
    C(1:sPreset.M,1:sPreset.M) = eye(sPreset.M)/sqrt(sPreset.n);
else
    b_maskDataTermCMatrix = false;
    C = RoMix(Phi, sPreset.gamma1, sPreset.gamma2, lambdaPhi, Ln, V, b_maskDataTermCMatrix);
end
VInt = (1/sqrt(interpRatio))*PhiInt*C;
VIntToCompare = VInt;
VRecPhi = Phi*C;
% ----------------------------------------------------------------------------------------------
% Interpolate with Nystrom
% ----------------------------------------------------------------------------------------------
if sPreset.b_flipSign
    VNys = FlipSign(VInt, VNys, sPreset.b_pairwiseFlipSign);
end
VNysToCompare = VNys;
% ----------------------------------------------------------------------------------------------
% Interpolate with Representer theorem
% ----------------------------------------------------------------------------------------------
b_normalizeAlpha = true;
mAlpha = LapRLS(W, V, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);

if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    CRep = diag(lambdaPhi)*Phi.'*mAlpha;
    PlotCoeffsMatrix(C, '${\bf C}$', CRep, ...
        '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
end
VRep = WTrainInt.'*mAlpha;
if sPreset.b_flipSign
    VRep = FlipSign(VInt, VRep, sPreset.b_pairwiseFlipSign);
end
VRepToCompare = VRep;
% ----------------------------------------------------------------------------------------------
% Calculate reference
% ----------------------------------------------------------------------------------------------
if sPreset.b_flipSign
    VRef = FlipSign(VInt, VRef, sPreset.b_pairwiseFlipSign);
    V = FlipSign(VRecPhi, V, sPreset.b_pairwiseFlipSign);
end
VRefToCompare = VRef;
VToCompare = V;
% ----------------------------------------------------------------------------------------------
% Plots
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    Cint = RoMix(PhiInt, sPreset.gamma1, sPreset.gamma2, lambdaPhi, LnRef, VRefToCompare, b_maskDataTermCMatrix);
    PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf RoMix}}$');
end
if sPlotParams.b_globalPlotEnable && sPreset.dim <= 3
    if sPreset.dim == 1
        plotInd = 0:4;
%         PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(end), ...
%             VRefToCompare, [], [], [], ['Reference vs. Representer theorem (N = ', num2str(N), ')'], ...
%             'VRep', 'v^{{\bf gt}}', VRepToCompare, 'v^{{\bf rep}}');
%         PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(end), ...
%             VRefToCompare, [], [], [], ['Reference vs. Nystrom (N = ', num2str(N), ')'], ...
%             'VNys', 'v^{{\bf gt}}', VNysToCompare, 'v^{{\bf nys}}');
        PlotEigenfuncvecScatter(sPlotParams, sPreset.verticesPDF, xInt, [], plotInd(1), plotInd(end), ...
            VRefToCompare, [], [], [], ['Eigenvectors interpolation (N = ', num2str(sPreset.N), ')'], ...
            'VInt', '\tilde{\psi}', VIntToCompare, '\tilde{\psi}^{{\bf RoMix}}');
    else
        plotInd = 1:4;
        [cData{1:numel(plotInd)*2}] = deal(xInt);
        cSigStr = [RepLegend('\\tilde{\\psi}', plotInd), RepLegend('\\tilde{\\psi}^{{\\bf RoMix}}', plotInd)];
        [cNumCircles{1:numel(plotInd)*2}] = deal(sPreset.N);
        [cMarkers{1:numel(plotInd)*2}] = deal('.');
        PlotGraphSignals(sPlotParams, ['Eigenvectors interpolation (N = ', num2str(sPreset.N), ')'], ...
            [sPreset.matrixForEigs, '_Eigs_',num2str(plotInd(1)), '_to_', num2str(plotInd(end)) ], cData, ...
            [mat2cell(VRefToCompare(:,plotInd+1),sPreset.N,ones(1,numel(plotInd))), mat2cell(VIntToCompare(:,plotInd+1),sPreset.N,ones(1,numel(plotInd)))], ...
            cSigStr, cNumCircles, cMarkers, [], [], [min(min(VRefToCompare(:,plotInd+1))), max(max(VRefToCompare(:,plotInd+1)))]);
%         cmap = PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(end), ...
%             VRefToCompare, [], [], [], ...
%             ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf gt}}');
%         PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(end), ...
%             VRepToCompare, [], [], [], ...
%             ['Representer theorem (N = ', num2str(N), ')'], 'VRep', 'v^{{\bf rep}}');
%         PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(end), ...
%             VNysToCompare, [], [], [], ...
%             ['Nystrom (N = ', num2str(N), ')'], 'VNys', 'v^{{\bf nys}}');
%         PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(end), ...
%             VIntToCompare, [], [], [], ...
%             ['Ours (N = ', num2str(N), ')'], 'VInt', 'v^{{\bf RoMix}}');
    end
end

% ----------------------------------------------------------------------------------------------
% Plot inner product of interpolated eigenvectors
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotInnerProductMatrices
    pltTitle = 'V - ${\bf V}^T {\bf V}$';
    figName = 'V';
    PlotInnerProductMatrix([], V, [], pltTitle, figName);
    pltTitle = 'VRef - ${\bf V}_{{\bf gt}}^T {\bf V}_{{\bf gt}}$';
    figName = 'VRef';
    PlotInnerProductMatrix([], VRef, [], pltTitle, figName);
    pltTitle = 'VInt - ${\bf V}_{{\bf RoMix}}^T {\bf V}_{{\bf RoMix}}$';
    figName = 'Vint';
    PlotInnerProductMatrix([], VInt, [], pltTitle, figName);
    pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
    figName = 'VNys';
    PlotInnerProductMatrix([], VNys, [], pltTitle, figName);
    pltTitle = 'VRep - ${\bf V}_{{\bf rep}}^T {\bf V}_{{\bf rep}}$';
    figName = 'VRep';
    PlotInnerProductMatrix([], VRep, [], pltTitle, figName);
%     pltTitle = '$\int \phi_i(x) \phi_j(x) p(x) dx = \Phi^T$diag(Pr)$\Phi$';
%     figName = 'Phi';
%     vPrTilde = sDistParams.vPr;
%     PlotInnerProductMatrix([], Phi, vPrTilde, pltTitle, figName);

end
end
