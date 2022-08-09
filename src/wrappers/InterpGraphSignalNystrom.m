function [mSigCnvrtRecV, mSigCnvrtNys, tTrain, tInt] = InterpGraphSignalNystrom(sPlotParams, sPreset, sDataset)

assert(~isempty(sDataset.sData.yt));
xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end
vLabeledInd = GetUnlabeledNodesMask(mSig);
% ------------------------------------------------------------------------------------------
% Build adjacency and perform eigs
% ------------------------------------------------------------------------------------------
[W, tTrainVec(1)] = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
ts = tic;
[V, Lambda] = eigs(W,sPreset.M);
adjLambda = diag(Lambda);
tTrainVec(2) = toc(ts);
% ------------------------------------------------------------------------------------------
% Reconstruction
% ------------------------------------------------------------------------------------------
ts = tic;
mSigCoeffsV = pinv(V(vLabeledInd,:))*mSig(vLabeledInd,:);
mSigRecV = V*mSigCoeffsV;
mSigCnvrtRecV = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecV);
tTrainVec(3) = toc(ts);
tTrain = sum(tTrainVec);

% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
ts = tic;
distTrainInt = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
WTrainInt = exp(-distTrainInt.^2/(2*sPreset.omega^2));
interpRatio = sPreset.N/sPreset.n;
lambdaNys = adjLambda*sqrt(interpRatio);
VNys = WTrainInt.'*V*diag(1./lambdaNys);
mSigNys = VNys*mSigCoeffsV;
mSigCnvrtNys = ConvertSignalByDataset(sPreset.verticesPDF, mSigNys);
tInt = toc(ts);

fprintf('Nystrom: train took %.2f sec, interp took %.2f sec\n', tTrain, tInt);

if sPlotParams.b_globalPlotEnable && sPreset.b_runGraphSignals
    mSigRef    = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigRefInt = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    PlotGraphSignalsWrapper([], sPreset, [], sDataset, mSigRef, mSigRefInt, mSigCnvrtRecV, mSigCnvrtNys, [], 'Nystrom')
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotInnerProductMatrices
    PlotInnerProductMatrix([], VNys, [], '${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$', 'VNys');
end

end