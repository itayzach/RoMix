function [VNysToCompare] = InterpEigenvectorsNystrom(sPlotParams, sPreset, sDataset, V, adjLambda, VInt)
xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;
% ----------------------------------------------------------------------------------------------
% Interpolate with Nystrom
% ----------------------------------------------------------------------------------------------
distLUBlockRUBlock = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
WTrainInt = exp(-distLUBlockRUBlock.^2/(2*sPreset.omega^2));
interpRatio = sPreset.N/sPreset.n;
lambdaNys = adjLambda*sqrt(interpRatio);
VNys = WTrainInt.'*V*diag(1./lambdaNys);
if sPreset.b_flipSign
    VNys = FlipSign(VInt, VNys, sPreset.b_pairwiseFlipSign);
end
VNysToCompare = VNys;

% ----------------------------------------------------------------------------------------------
% Plot inner product of interpolated eigenvectors
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotInnerProductMatrices
    pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
    figName = 'VNys';
    PlotInnerProductMatrix([], VNys, [], pltTitle, figName);
   
end
end
