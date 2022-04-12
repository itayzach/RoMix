function mSignalOut = ConvertSignalByDataset(verticesPDF, mSignalIn)
if ismember(verticesPDF, {'TwoMoons', 'TwoSpirals'})
    mSignalOut = sign(mSignalIn);
elseif ismember(verticesPDF, {'USPS', 'MNIST', 'MnistLatentVAE'})  
%     mSig       = sigmoid(mSig);
%     mSigRecPhi = lowpass(mSigRecPhi,0.01);
    [~, mSignalOut] = max(mSignalIn,[],2);
else
    mSignalOut = mSignalIn;
end
end