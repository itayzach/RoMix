function [mSigCnvrtRecV, mSigCnvrtNys] = InterpGraphSignalNystrom(sPreset, sDataset, V, adjLambda)

assert(~isempty(sDataset.sData.yt));
xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;

% ------------------------------------------------------------------------------------------
% Reconstruction
% ------------------------------------------------------------------------------------------
mSig = sDataset.sData.y;
mSigCoeffsV = V'*mSig; % same as pinv(V)*sig...

% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
distLUBlockRUBlock = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
WTrainInt = exp(-distLUBlockRUBlock.^2/(2*sPreset.omega^2));
interpRatio = sPreset.N/sPreset.n;
lambdaNys = adjLambda*sqrt(interpRatio);
VNys = WTrainInt.'*V*diag(1./lambdaNys);

mSigRecV   = V*mSigCoeffsV;
mSigNys    = VNys*mSigCoeffsV;

mSigCnvrtRecV   = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecV);
mSigCnvrtNys    = ConvertSignalByDataset(sPreset.verticesPDF, mSigNys);

end