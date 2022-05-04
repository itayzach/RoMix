function [mSigCnvrtRecKnn, mSigCnvrtKnn, vTrainTime, vIntTime] = InterpGraphSignalKnn(sPlotParams, sPreset, sDataset)

assert(~isempty(sDataset.sData.yt));
xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end

% ------------------------------------------------------------------------------------------
% Reconstruction
% ------------------------------------------------------------------------------------------
%[idx, D] = knnsearch(X1, X2, 'K', k);
%WNN = NearestNeighborsAdjacency(xTrain, xTrain, sPreset.knn, sPreset.nnValue);
ts = tic;
mSigKnn = mSig;
J = GetUnlabeledNodesMask(mSig);
mSigCnvrtRecKnn = ConvertSignalByDataset(sPreset.verticesPDF, mSigKnn);
vTrainTime = toc(ts);

% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
ts = tic;
n = size(xTrain,1);
N = size(xInt,1);
% idx = knnsearch(xTrain, xInt(n+1:end,:), 'K', sPreset.knn);
% mSigKnnInt(1:n,:) = mSig;
% for i = 1:N-n
%     mSigKnnInt(i+n,:) = (1/sPreset.knn)*sum(mSig(idx(i,:)));
% end
mSigKnnInt(1:n,:) = mSig;
for i = 1:N-n
    idx = knnsearch(xTrain, xInt(i+n,:), 'K', sPreset.knn);
    mSigKnnInt(i+n,:) = (1/sPreset.knn)*sum(mSig(idx,:));
end

mSigCnvrtKnn = ConvertSignalByDataset(sPreset.verticesPDF, mSigKnnInt);
vIntTime = toc(ts);
fprintf('KNN: train took %.2f sec, interp took %.2f sec\n', vTrainTime, vIntTime);

if sPlotParams.b_globalPlotEnable && sPreset.b_runGraphSignals
    mSigRef    = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    mSigRefInt = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
    PlotGraphSignalsWrapper([], sPreset, [], sDataset, mSigRef, mSigRefInt, mSigCnvrtRecKnn, mSigCnvrtKnn, [], 'kNN')
end

end