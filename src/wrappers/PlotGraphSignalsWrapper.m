function PlotGraphSignalsWrapper(sPlotParams, sPreset, sKernelParams, sDataset, mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt, mSigCoeffsPhi, methodName)

xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
if isfield(sDataset.sData, 'ymasked')
    mSig = sDataset.sData.ymasked;
else
    mSig = sDataset.sData.y;
end
vLabeledInd = GetUnlabeledNodesMask(mSig);
% vUnlabeledInd = find(~vLabeledFlag);
% xTrainLabeled = sDataset.sData.x(vLabeledInd,:);
% xTrainUnlabeled = sDataset.sData.x(vUnlabeledInd,:);
% nUnlabeled = size(xTrainUnlabeled,1);
% nLabeled = size(xTrainLabeled,1);
% mSigLabeled = mSig(vLabeledInd,:);

% ------------------------------------------------------------------------------------------
% Plots
% ------------------------------------------------------------------------------------------
sigIndToPlot = min(7, size(mSigCnvrtRecRef,2));

vSigRecRef = mSigCnvrtRecRef(:,sigIndToPlot);
vSigIntRef = mSigCnvrtIntRef(:,sigIndToPlot);
vSigRec    = mSigCnvrtRec(:,sigIndToPlot);
vSigInt    = mSigCnvrtInt(:,sigIndToPlot);

if ismember(sPreset.verticesPDF, {'BulgariBeacons'})
    grid.lat = sDataset.sData.vGridLat;
    grid.lon = sDataset.sData.vGridLon;
    PlotGraphSignals(sPlotParams, ['Graph signal reconstruction \& interpolation'], 'Interpolation', ...
       {sDataset.sData.mGrid, sDataset.sData.mGrid}, ...
       {vSigRecRef, vSigRec}, ... 
       {'$s$', ['$s^{{\bf ' methodName '}}$']}, ...
       {1:size(xTrain,1), 1:size(xTrain,1)}, [], [], [], [], grid, sDataset.sData.tx, sDataset.sData.sense);

    if size(sDataset.sData.x,2) <= 3
        PlotGraphSignals(sPlotParams, ['Graph signal reconstruction \& interpolation'], 'Interpolation', ...
            {sDataset.sData.x, sDataset.sData.xt, sDataset.sData.x, sDataset.sData.xt, sDataset.sData.x}, ...
            {vSigRecRef, vSigIntRef, vSigRec, vSigInt, sDataset.sData.ymasked }, ... cellfun(@db,{vSigRecRef, vSigIntRef, vSigRec, vSigInt}, 'UniformOutput',false), ...
            {'$s$', '$\tilde{s}$', ['$s^{{\bf ' methodName '}}$'], ['$\tilde{s}^{{\bf ' methodName '}}$'], 'ymasked'}, ...
            {1:size(xTrain,1), 1:size(xTrain,1), 1:size(xTrain,1), 1:size(xTrain,1), 1:size(xTrain,1)});
    end
end

if sPreset.dim <=3
    if sPreset.dim == 1
        colorOrder = get(gca, 'ColorOrder');
        PlotGraphSignals(sPlotParams, ['Graph signal reconstruction'], 'Reconstruction', ...
            {xTrain, xTrain}, ...
            {vSigRecRef, vSigRec}, ...
            {'$s$', ['$s^{{\bf ' methodName '}}$']}, ...
            {vLabeledInd, vLabeledInd}, {'o', '.'}, mat2cell(colorOrder(1:2,:),[1 1], 3));
        PlotGraphSignals(sPlotParams, ['Graph signal interpolation'], 'Interpolation', ...
            {xInt, xInt}, ...
            {vSigIntRef, vSigInt}, ...
            {'$\tilde{s}$', ['$\tilde{s}^{{\bf ' methodName '}}$']}, ...
            {vLabeledInd, vLabeledInd}, {'o', '.'}, mat2cell(colorOrder(3:4,:),[1 1], 3));
    else
        PlotGraphSignals(sPlotParams, ['Graph signal reconstruction \& interpolation'], 'Interpolation', ...
            {xTrain, xInt, xTrain, xInt}, ...
            {vSigRecRef, vSigIntRef, vSigRec, vSigInt}, ...
            {'$s$', '$\tilde{s}$', ['$s^{{\bf ' methodName '}}$'], ['$\tilde{s}^{{\bf ' methodName '}}$']}, ...
            {vLabeledInd, vLabeledInd, vLabeledInd, vLabeledInd});

    end
    if ~isempty(sPlotParams) && sPlotParams.b_plotGmmSignal
        % GMM
        nGmmPoints = 1000;
        [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
        [PhiGmm, ~] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xGmm);
        mSigGmm = PhiGmm*mSigCoeffsPhi;
        vSigGmm = mSigGmm(:,sigIndToPlot);
        if sPreset.dim == 2
            xylim(1) = min(xTrain(:,1));
            xylim(2) = max(xTrain(:,1));
            xylim(3) = min(xTrain(:,2));
            xylim(4) = max(xTrain(:,2));
        else
            xylim = [];
        end
        PlotGraphSignals(sPlotParams, 'GMM Graph signal', 'GMM', {xGmm}, {vSigGmm}, ...
            {'$s^{{\bf gmm}}$'}, {nGmmPoints}, {'.'}, xylim);
    end
end

if ismember(sPreset.verticesPDF, {'TwoMoons'})
    % make sure mSigCnvrtRecPhi instead of mSigRecPhi
    warning('TODO: Update PlotTwoMoonsRoMix to support other methods')
    if(strcmp(methodName,'RoMix'))
        PlotTwoMoonsRoMix(sPlotParams, sDataset, sKernelParams, mSigCnvrtRec, mSigCoeffsPhi, sPreset.gamma1, sPreset.gamma2);
    end

elseif ismember(sPreset.verticesPDF, {'USPS', 'MNIST'}) && ~isempty(sKernelParams)
    b_transpose = strcmp(sPreset.verticesPDF, 'MNIST') && sPreset.dim == 28*28;
    b_minusOne = strcmp(sPreset.verticesPDF, 'MNIST');
    nDigitsToPlot = 50;    
    assert(sPreset.n + nDigitsToPlot < sPreset.N)
    vSamples = round(linspace(sPreset.n+1,sPreset.N,nDigitsToPlot));
    if isfield(sKernelParams.sDistParams, 'vae')
        % Transform from z to x
        xDigits = reshape(double(pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=xInt,vae=sKernelParams.sDistParams.vae)),[],28*28);
    else
        xDigits = xInt;
    end
    
    figTitle = ['Prediction on unseen $', num2str(length(vSamples)), '$ images'];
    figName = 'interpolation_signals';
    PlotDigits(sPlotParams, xDigits(vSamples,:), mSigCnvrtInt(vSamples)-b_minusOne, b_transpose, figTitle, figName);

    nBadPredict = 50;
    vBadPredictInd = find(mSigCnvrtInt ~= mSigCnvrtIntRef, nBadPredict);
    nBadPredict = min(length(vBadPredictInd));
    figTitle = ['Wrong prediction on given+unseen $n = ', num2str(nBadPredict), '$ points'];
    PlotDigits([], xDigits(vBadPredictInd,:), mSigCnvrtInt(vBadPredictInd)-b_minusOne, b_transpose, figTitle, 'Wrong predictions');
    
    nGmmPoints = 50;
    [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
    if isfield(sKernelParams.sDistParams, 'vae')
        % Transform from z to x
        xDigits = reshape(double(pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=xGmm,vae=sKernelParams.sDistParams.vae)),[],28*28);
    else
        xDigits = xGmm;
    end
    PhiGmm = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xGmm);
    mSigGmm = PhiGmm*mSigCoeffsPhi;
    mSigCnvrtGmm = ConvertSignalByDataset(sPreset.verticesPDF, mSigGmm);
    figTitle = [num2str(nGmmPoints), ' generated points'];
    figName = 'generated_signals';
    PlotDigits(sPlotParams, xDigits, mSigCnvrtGmm-b_minusOne, b_transpose, figTitle, figName);

end
end
