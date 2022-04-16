function PlotGraphSignalsWrapper(sPlotParams, sPreset, sKernelParams, sDataset, mSigCnvrt, mSigCnvrtRef, mSigCnvrtRecPhi, mSigCnvrtInt, mSigCoeffsPhi)

xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;

% ------------------------------------------------------------------------------------------
% Plots
% ------------------------------------------------------------------------------------------
sigIndToPlot     = min(7, size(mSigCnvrt,2));

vSig       = mSigCnvrt(:,sigIndToPlot);
vSigRef    = mSigCnvrtRef(:,sigIndToPlot);
vSigRecPhi = mSigCnvrtRecPhi(:,sigIndToPlot);
vSigInt    = mSigCnvrtInt(:,sigIndToPlot);

if sPreset.dim <=3
    if sPreset.dim == 1
        colorOrder = get(gca, 'ColorOrder');
        PlotGraphSignals(sPlotParams, ['Graph signal reconstruction'], 'Reconstruction', ...
            {xTrain, xTrain}, ...
            {vSig, vSigRecPhi}, ...
            {'$s$', '$s^{{\bf RoMix}}$'}, ...
            {sPreset.n, sPreset.n}, {'o', '.'}, mat2cell(colorOrder(1:2,:),[1 1], 3));
        PlotGraphSignals(sPlotParams, ['Graph signal interpolation'], 'Interpolation', ...
            {xInt, xInt}, ...
            {vSigRef, vSigInt}, ...
            {'$\tilde{s}$', '$\tilde{s}^{{\bf RoMix}}$'}, ...
            {sPreset.n, sPreset.n}, {'o', '.'}, mat2cell(colorOrder(3:4,:),[1 1], 3));
    else
        PlotGraphSignals(sPlotParams, ['Graph signal reconstruction \& interpolation'], 'Interpolation', ...
            {xTrain, xInt, xTrain, xInt}, ...
            {vSig, vSigRef, vSigRecPhi, vSigInt}, ...
            {'$s$', '$\tilde{s}$', '$s^{{\bf RoMix}}$', '$\tilde{s}^{{\bf RoMix}}$'}, ...
            {sPreset.n, sPreset.n, sPreset.n, sPreset.n});
    end
    if sPlotParams.b_plotGmmSignal
        % GMM
        nGmmPoints = 1000;
        [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
        [PhiGmm, ~] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xGmm, sPreset.b_normalizePhi);
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
    PlotTwoMoonsRoMix(sPlotParams, sDataset, sKernelParams, mSigCnvrtRecPhi, mSigCoeffsPhi, sPreset.gamma1, sPreset.gamma2, sPreset.b_normalizePhi);

elseif ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})
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

    vBadPredictInd = find(mSigCnvrtInt ~= mSigCnvrtRef);
    nBadPredict = min(50,length(vBadPredictInd));
    figTitle = ['Wrong prediction on given+unseen $n = ', num2str(nBadPredict), '$ points'];
    PlotDigits([], xDigits(vBadPredictInd(1:nBadPredict),:), mSigCnvrtInt(vBadPredictInd(1:nBadPredict))-b_minusOne, b_transpose, figTitle);
    
    nGmmPoints = 50;
    [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
    if isfield(sKernelParams.sDistParams, 'vae')
        % Transform from z to x
        xDigits = reshape(double(pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=xGmm,vae=sKernelParams.sDistParams.vae)),[],28*28);
    else
        xDigits = xGmm;
    end
    PhiGmm = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xGmm, sPreset.b_normalizePhi);
    mSigGmm = PhiGmm*mSigCoeffsPhi;
    mSigCnvrtGmm = ConvertSignalByDataset(sPreset.verticesPDF, mSigGmm);
    figTitle = [num2str(nGmmPoints), ' generated points'];
    figName = 'generated_signals';
    PlotDigits(sPlotParams, xDigits, mSigCnvrtGmm-b_minusOne, b_transpose, figTitle, figName);

end
end