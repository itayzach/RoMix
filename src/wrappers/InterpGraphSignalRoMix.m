function [mSigCnvrtRecPhi, mSigCnvrt, mSigCnvrtInt, mSigCnvrtRef, C] = ...
    InterpGraphSignalRoMix(sPreset, sDataset, Phi, lambdaPhi, PhiInt, Ln)
assert(~isempty(sDataset.sData.yt));
mSig = sDataset.sData.y;
mSigRef = sDataset.sData.yt;
% ------------------------------------------------------------------------------------------
% Coeffs
% ------------------------------------------------------------------------------------------
invLambda = diag(1./lambdaPhi);
if isfield(sDataset.sData, 'ymasked')
    mSigMasked = sDataset.sData.ymasked;
    C = RoMix(Phi, sPreset.gamma1, sPreset.gamma2, invLambda, Ln, mSigMasked, sPreset.b_maskDataFitTerm);
else
    C = RoMix(Phi, sPreset.gamma1, sPreset.gamma2, invLambda, Ln, mSig, sPreset.b_maskDataFitTerm);
end

% ------------------------------------------------------------------------------------------
% Signals
% ------------------------------------------------------------------------------------------
mSigRecPhi = Phi*C;
mSigInt    = PhiInt*C;

mSigCnvrt       = ConvertSignalByDataset(sPreset.verticesPDF, mSig);
mSigCnvrtRef    = ConvertSignalByDataset(sPreset.verticesPDF, mSigRef);
mSigCnvrtRecPhi = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecPhi);
mSigCnvrtInt    = ConvertSignalByDataset(sPreset.verticesPDF, mSigInt);

end