function [VIntToCompare, VNysToCompare, VRepToCompare, VRefToCompare] = ...
    InterpEigenvectors(sPlotParams, sPreset, sDataset, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, VRef, W, WRef, D, DRef, Ln, LnRef)
dim                = sPreset.dim;
n                  = sPreset.n;
N                  = sPreset.N;
verticesPDF        = sPreset.verticesPDF;
M                  = sPreset.M;
MTilde             = sPreset.MTilde;
gamma1             = sPreset.gamma1;
gamma2             = sPreset.gamma2;
gamma1Rep          = sPreset.gamma1Rep;
gamma2Rep          = sPreset.gamma2Rep;
b_forceCtoIdentity = sPreset.b_forceCtoIdentity;
b_flipSign         = sPreset.b_flipSign;
b_pairwiseFlipSign = sPreset.b_pairwiseFlipSign;
interpRatio        = N/n;
if dim == 1
    plotInd = [0,min(4,M-1)];
else
    plotInd = [0,min(8,M-1)];
end
xInt = sDataset.sData.xt;
% ----------------------------------------------------------------------------------------------
% Interpolate with our method
% ----------------------------------------------------------------------------------------------
if b_forceCtoIdentity
    C = zeros(MTilde, M);
    C(1:M,1:M) = eye(M)/sqrt(n);
else
    b_maskDataTermCMatrix = false;
    invLambda = diag(1./lambdaPhi);
    C = EigsRLS(Phi, gamma1, gamma2, invLambda, Ln, V, b_maskDataTermCMatrix);
end
VInt = (1/sqrt(interpRatio))*PhiInt*C;
VIntToCompare = VInt;
% ----------------------------------------------------------------------------------------------
% Interpolate with Nystrom
% ----------------------------------------------------------------------------------------------
if b_flipSign
    VNys = FlipSign(VInt, VNys, b_pairwiseFlipSign);
end
VNysToCompare = VNys;
% ----------------------------------------------------------------------------------------------
% Interpolate with Representer theorem
% ----------------------------------------------------------------------------------------------
b_normalizeAlpha = true;
mAlpha = LapRLS(W, V, Ln, gamma1Rep, gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);

if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    CRep = diag(lambdaPhi)*Phi.'*mAlpha;
    PlotCoeffsMatrix(C, '${\bf C}$', CRep, ...
        '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
end
VRep = WTrainInt.'*mAlpha;
if b_flipSign
    VRep = FlipSign(VInt, VRep, b_pairwiseFlipSign);
end
VRepToCompare = VRep;
% ----------------------------------------------------------------------------------------------
% Calculate reference
% ----------------------------------------------------------------------------------------------
if b_flipSign
    VRef = FlipSign(VInt, VRef, b_pairwiseFlipSign);
end
VRefToCompare = VRef;

% ----------------------------------------------------------------------------------------------
% Plots
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    Cint = EigsRLS(PhiInt, gamma1, gamma2, invLambda, LnRef, VRefToCompare, b_maskDataTermCMatrix);
    PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf int}}$');
end
if sPlotParams.b_globalPlotEnable && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
    if dim == 1
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VRefToCompare, [], [], [], ['Reference vs. Representer theorem (N = ', num2str(N), ')'], ...
            'VRep', 'v^{{\bf ref}}', VRepToCompare, 'v^{{\bf rep}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VRefToCompare, [], [], [], ['Reference vs. Nystrom (N = ', num2str(N), ')'], ...
            'VNys', 'v^{{\bf ref}}', VNysToCompare, 'v^{{\bf nys}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VRefToCompare, [], [], [], ['Reference vs. Ours (N = ', num2str(N), ')'], ...
            'VInt', 'v^{{\bf ref}}', VIntToCompare, 'v^{{\bf int}}');
    else
        cmap = PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VRefToCompare, [], [], [], ...
            ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf ref}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
            VRepToCompare, [], [], [], ...
            ['Representer theorem (N = ', num2str(N), ')'], 'VRep', 'v^{{\bf rep}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
            VNysToCompare, [], [], [], ...
            ['Nystrom (N = ', num2str(N), ')'], 'VNys', 'v^{{\bf nys}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
            VIntToCompare, [], [], [], ...
            ['Ours (N = ', num2str(N), ')'], 'VInt', 'v^{{\bf int}}');
    end
end

% ----------------------------------------------------------------------------------------------
% Plot inner product of interpolated eigenvectors
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotInnerProductMatrices
    pltTitle = 'VRef - ${\bf V}_{{\bf ref}}^T {\bf V}_{{\bf ref}}$';
    figName = 'VRef';
    PlotInnerProductMatrix([], VRef, [], pltTitle, figName);
    pltTitle = 'VInt - ${\bf V}_{{\bf int}}^T {\bf V}_{{\bf int}}$';
    figName = 'Vint';
    PlotInnerProductMatrix([], VInt, [], pltTitle, figName);
    pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
    figName = 'VNys';
    PlotInnerProductMatrix([], VNys, [], pltTitle, figName);
    pltTitle = 'VRep - ${\bf V}_{{\bf rep}}^T {\bf V}_{{\bf rep}}$';
    figName = 'VRep';
    PlotInnerProductMatrix([], VRep, [], pltTitle, figName);
    %             pltTitle = '$\int \phi_i(x) \phi_j(x) p(x) dx = \Phi^T$diag(Pr)$\Phi$';
    %             figName = 'Phi';
    %             vPrTilde = sDistParams.vPr;
    %             PlotInnerProductMatrix([], Phi, vPrTilde, pltTitle, figName);

end
end