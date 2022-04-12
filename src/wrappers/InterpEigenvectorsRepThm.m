function [VRepToCompare] = InterpEigenvectorsRepThm(sPlotParams, sPreset, sDataset, W, V, VInt, Ln)
interpRatio = sPreset.N/sPreset.n;

xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;

% ----------------------------------------------------------------------------------------------
% Interpolate with Representer theorem
% ----------------------------------------------------------------------------------------------
distLUBlockRUBlock = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
WTrainInt = exp(-distLUBlockRUBlock.^2/(2*sPreset.omega^2));

b_normalizeAlpha = true;
mAlpha = LapRLS(W, V, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);

VRep = WTrainInt.'*mAlpha;
if sPreset.b_flipSign
    VRep = FlipSign(VInt, VRep, sPreset.b_pairwiseFlipSign);
end
VRepToCompare = VRep;

% ----------------------------------------------------------------------------------------------
% Plot inner product of interpolated eigenvectors
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotInnerProductMatrices
    pltTitle = 'VRep - ${\bf V}_{{\bf rep}}^T {\bf V}_{{\bf rep}}$';
    figName = 'VRep';
    PlotInnerProductMatrix([], VRep, [], pltTitle, figName);

end
end
