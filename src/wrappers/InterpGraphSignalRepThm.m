function [mSigCnvrtRecRep, mSigCnvrtRep] = InterpGraphSignalRepThm(sPreset, sDataset, W, Ln)

interpRatio = sPreset.N/sPreset.n;
assert(~isempty(sDataset.sData.yt));
b_normalizeAlpha = false;

xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;

% ------------------------------------------------------------------------------------------
% Coeffs
% ------------------------------------------------------------------------------------------
if isfield(sDataset.sData, 'ymasked')
    mSigMasked = sDataset.sData.ymasked;
    mAlpha = LapRLS(W, mSigMasked, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
else
    mSig = sDataset.sData.y;
    mAlpha = LapRLS(W, mSig, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
end

% ------------------------------------------------------------------------------------------
% Interpolation
% ------------------------------------------------------------------------------------------
distLUBlockRUBlock = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
WTrainInt = exp(-distLUBlockRUBlock.^2/(2*sPreset.omega^2));
mSigRecRep = W.'*mAlpha;
mSigRep    = WTrainInt.'*mAlpha;

mSigCnvrtRecRep = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecRep);
mSigCnvrtRep    = ConvertSignalByDataset(sPreset.verticesPDF, mSigRep);

end