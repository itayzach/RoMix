function [mSigCnvrtRecPw, mSigCnvrtPw, trainTime, intTime] = InterpGraphSignalPesenson(sPlotParams, sPreset, sDataset)

xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
n = size(xTrain,1);
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end
vLabeledFlag = diag(GetUnlabeledNodesMask(mSig));
vLabeledInd = find(vLabeledFlag);
mSigLabeled = mSig(vLabeledInd,:);
% ------------------------------------------------------------------------------------------
% Reconstruction (on unlabeled training nodes)
% ------------------------------------------------------------------------------------------
[W, tTrainVec(1)] = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);

ts = tic;
G = gsp_graph(W, xTrain);
G.lap_type='normalized';
G = gsp_create_laplacian(G);
G = gsp_estimate_lmax(G);
mSigPw = gsp_interpolate(G, mSigLabeled, vLabeledInd, sPreset.sPwParams);
mSigCnvrtRecPw = ConvertSignalByDataset(sPreset.verticesPDF, mSigPw);
tTrainVec(2) = toc(ts);
trainTime = sum(tTrainVec);
% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
[WRef, tIntVec(1)] = CalcAdjacency(xInt, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
ts = tic;
GInt = gsp_graph(WRef, xInt);
GInt.lap_type='normalized';
GInt = gsp_create_laplacian(GInt);
GInt = gsp_estimate_lmax(GInt);
mSigPwInt = gsp_interpolate(GInt, mSigLabeled, vLabeledInd, sPreset.sPwParams);
mSigCnvrtPw = ConvertSignalByDataset(sPreset.verticesPDF, mSigPwInt);
tIntVec(2) = toc(ts);
intTime = sum(tIntVec);

fprintf('PW: train took %.2f sec, interp took %.2f sec\n', trainTime, intTime);

if sPlotParams.b_globalPlotEnable && sPreset.b_runGraphSignals
    mSigRef    = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigRefInt = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    PlotGraphSignalsWrapper([], sPreset, [], sDataset, mSigRef, mSigRefInt, mSigCnvrtRecPw, mSigCnvrtPw, [], 'PW')
end
end