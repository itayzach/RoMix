function [mSigCnvrtRecRep, mSigCnvrtRep, tTrain, tInt, mAlpha] = InterpGraphSignalRepThm(sPlotParams, sPreset, sDataset)

b_normalizeAlpha = false;

xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end
% ------------------------------------------------------------------------------------------
% Build adjacency and find alpha
% ------------------------------------------------------------------------------------------
[W, tTrainVec(1), ~, ~, Ln, tTrainVec(2)] = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
[mAlpha, tTrainVec(3)] = LapRLS(W, mSig, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, sPreset.interpRatio, b_normalizeAlpha);
% ------------------------------------------------------------------------------------------
% Reconstruction
% ------------------------------------------------------------------------------------------
ts = tic;
mSigRecRep = W.'*mAlpha;
mSigCnvrtRecRep = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecRep);
tTrainVec(3) = toc(ts);
tTrain = sum(tTrainVec);
% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
ts = tic;
distTrainInt = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
WTrainInt = exp(-distTrainInt.^2/(2*sPreset.omega^2));
mSigRep = WTrainInt.'*mAlpha;
mSigCnvrtRep = ConvertSignalByDataset(sPreset.verticesPDF, mSigRep);
tInt = toc(ts);

fprintf('Rep. Thm.: train took %.2f sec, interp took %.2f sec\n', tTrain, tInt);

if sPlotParams.b_globalPlotEnable && sPreset.b_runGraphSignals
    mSigRef    = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigRefInt = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    PlotGraphSignalsWrapper(sPlotParams, sPreset, [], sDataset, mSigRef, mSigRefInt, mSigCnvrtRecRep, mSigCnvrtRep, [], [], [], 'Rep. Thm.')
end
end