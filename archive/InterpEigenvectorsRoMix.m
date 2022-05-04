function [VIntToCompare, VRecPhi, VToCompare, VRefToCompare, C] = ...
    InterpEigenvectorsRoMix(sPlotParams, sPreset, Phi, lambdaPhi, PhiInt, V, VRef, Ln)
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
% Calculate reference
% ----------------------------------------------------------------------------------------------
if sPreset.b_flipSign
    VRef = FlipSign(VInt, VRef, sPreset.b_pairwiseFlipSign);
    V = FlipSign(VRecPhi, V, sPreset.b_pairwiseFlipSign);
end
VRefToCompare = VRef;
VToCompare = V;

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
end
end
