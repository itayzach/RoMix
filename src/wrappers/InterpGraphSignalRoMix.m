function [mSigCnvrtRec, mSigCnvrtInt, tTrain, tInt, C] = InterpGraphSignalRoMix(sPlotParams, sPreset, sDataset)
xTrain = sDataset.sData.x; 
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end

% 1. Estimate distribution parameters (GMM)
[sDistParams, tTrainVec(1:2)] = EstimateDistributionParameters(xTrain, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);
if isfield(sDataset, 'vae'), sDistParams.vae = sDataset.vae; end
PlotRoMixDatasetAndGmm(sPlotParams, sPreset, sDataset, sDistParams);

% 2. Calculate Phi(xTrain) and lambdaPhi
[sKernelParams, tTrainVec(3)] = CalcKernelParams(sDistParams, sPreset.omegaTilde);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex, tTrainVec(4)] = CalcAnalyticEigenvalues(sPreset.MTilde, sKernelParams);
[Phi, lambdaPhi, tTrainVec(4)] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xTrain);

% 3. Find C
if sPreset.gamma2 > 0
    [~, ~, ~, ~, Ln, tTrainVec(5)] = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
else
    Ln = [];
    tTrainVec(5) = 0;
end
if isfield(sDataset.sData, 'vLabeledInd')
    % Make sure we are consistent...
    % Some of the presets use vLabeledInd, and some sDataset.sData.ymasked
    assert(isequal(GetUnlabeledNodesMask(mSig),sort(sDataset.sData.vLabeledInd)) ...
        && numel(GetUnlabeledNodesMask(mSig)) == sPreset.nLabeled);
end
[C, tTrainVec(6)] = RoMix(Phi, sPreset.gamma1, sPreset.gamma2, lambdaPhi, Ln, mSig);

% 4. Project (reconstruct)
mSigRecPhi = Phi*C;
mSigCnvrtRec = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecPhi);

% 5. Calculate Phi(xInt) and interpolate
[PhiInt, ~, tIntVec(1)] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xInt);
ts = tic;
mSigInt = PhiInt*C;
mSigCnvrtInt = ConvertSignalByDataset(sPreset.verticesPDF, mSigInt);
tIntVec(2) = toc(ts);

% 6. Take time
tTrain = sum(tTrainVec);
tInt = sum(tIntVec);
fprintf('RoMix: train took ');
fprintf('%.2f + ', tTrainVec(1:end-1))
fprintf('%.2f = %.2f sec, interp took %.2f sec\n', tTrainVec(end), tTrain, tInt);

% 7. Plot
PlotRoMixEigsAnalysis(sPlotParams, sPreset, sDataset, sDistParams, sKernelParams, C, Phi, lambdaPhi, PhiInt);
PlotRoMixGraphSignals(sPlotParams, sPreset, sDataset, sKernelParams, mSigCnvrtRec, mSigCnvrtInt, C);

end

