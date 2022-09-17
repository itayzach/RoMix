function [mSigCnvrtKnnRec, mSigCnvrtKnnInt, trainTime, intTime] = InterpGraphSignalKnn(sPlotParams, sPreset, sDataset)

assert(~isempty(sDataset.sData.yt));
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end
vLabeledInd = GetUnlabeledNodesMask(mSig);
vUnlabeledInd = setdiff(1:sPreset.n,vLabeledInd);
xTrainLabeled = sDataset.sData.x(vLabeledInd,:);
xTrainUnlabeled = sDataset.sData.x(vUnlabeledInd,:);
nUnlabeled = size(xTrainUnlabeled,1);
nLabeled = size(xTrainLabeled,1);
mSigLabeled = mSig(vLabeledInd,:);
% ------------------------------------------------------------------------------------------
% Reconstruction (on unlabeled training nodes)
% ------------------------------------------------------------------------------------------
ts = tic;
mSigKnnRec = mSig;
for i = 1:nUnlabeled
    [idx, D] = knnsearch(xTrainLabeled, xTrainUnlabeled(i,:), 'K', sPreset.knn);
    w = D'./sum(D);
    mSigKnnRec(vUnlabeledInd(i),:) = sum(w.*mSigLabeled(idx,:));
end
mSigCnvrtKnnRec = ConvertSignalByDataset(sPreset.verticesPDF, mSigKnnRec);
trainTime = toc(ts);

% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
ts = tic;
n = nLabeled+nUnlabeled;
assert(n == sPreset.n);
N = size(xInt,1);
mSigKnnInt(1:n,:) = mSigKnnRec;
for i = 1:N-n
    [idx, D] = knnsearch(xTrainLabeled, xInt(i+n,:), 'K', sPreset.knn);
    w = D'./sum(D);
    mSigKnnInt(i+n,:) = sum(w.*mSigLabeled(idx,:));
end
mSigCnvrtKnnInt = ConvertSignalByDataset(sPreset.verticesPDF, mSigKnnInt);
intTime = toc(ts) + trainTime;
fprintf('KNN: train took %.2f sec, interp took %.2f sec\n', trainTime, intTime);

if sPlotParams.b_globalPlotEnable && sPreset.b_runGraphSignals
    mSigRef    = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigRefInt = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    PlotGraphSignalsWrapper(sPlotParams, sPreset, [], sDataset, mSigRef, mSigRefInt, mSigCnvrtKnnRec, mSigCnvrtKnnInt, [], 'kNN')
end

end
