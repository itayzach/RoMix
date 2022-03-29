function [mSigCnvrtRecPhi, mSigCnvrtRecV, mSigCnvrtRecRep, mSigCnvrt, mSigCnvrtInt, mSigCnvrtNys, mSigCnvrtRep, mSigCnvrtRef] = ...
    InterpGraphSignal(sPlotParams, sPreset, sDataset, sKernelParams, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, W, WRef, D, DRef, Ln, LnRef)
interpRatio        = sPreset.N/sPreset.n;

assert(~isempty(sDataset.sData.yt));
xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
mSig = sDataset.sData.y;
mSigRef = sDataset.sData.yt;
% ------------------------------------------------------------------------------------------
% Coeffs
% ------------------------------------------------------------------------------------------
mSigCoeffsV = V'*mSig; % same as pinv(V)*sig...

invLambda = diag(1./lambdaPhi);
if isfield(sDataset.sData, 'ymasked')
    mSigMasked = sDataset.sData.ymasked;
    mSigCoeffsPhi = EigsRLS(Phi, sPreset.gamma1, sPreset.gamma2, invLambda, Ln, mSigMasked, sPreset.b_maskDataFitTerm);
    b_normalizeAlpha = false;
    mAlpha = LapRLS(W, mSigMasked, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
else
    mSigCoeffsPhi = EigsRLS(Phi, sPreset.gamma1, sPreset.gamma2, invLambda, Ln, mSig, sPreset.b_maskDataFitTerm);
    b_normalizeAlpha = false;
    mAlpha = LapRLS(W, mSig, Ln, sPreset.gamma1Rep, sPreset.gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
end

% Just for reference
mSigRefCoeffsPhi = EigsRLS(PhiInt, sPreset.gamma1, sPreset.gamma2, invLambda, LnRef, mSigRef, sPreset.b_maskDataFitTerm);
% ------------------------------------------------------------------------------------------
% Signals
% ------------------------------------------------------------------------------------------
mSigRecV   = V*mSigCoeffsV;
mSigRecPhi = Phi*mSigCoeffsPhi;
mSigRecRep = W.'*mAlpha;
mSigInt    = PhiInt*mSigCoeffsPhi;
mSigNys    = VNys*mSigCoeffsV;
mSigRep    = WTrainInt.'*mAlpha;
mSigIntRef = PhiInt*mSigRefCoeffsPhi;

mSigCnvrt       = ConvertSignalByDataset(sPreset.verticesPDF, mSig);
mSigCnvrtRef    = ConvertSignalByDataset(sPreset.verticesPDF, mSigRef);
mSigCnvrtRecV   = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecV);
mSigCnvrtRecPhi = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecPhi);
mSigCnvrtRecRep = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecRep);
mSigCnvrtInt    = ConvertSignalByDataset(sPreset.verticesPDF, mSigInt);
mSigCnvrtNys    = ConvertSignalByDataset(sPreset.verticesPDF, mSigNys);
mSigCnvrtRep    = ConvertSignalByDataset(sPreset.verticesPDF, mSigRep);
mSigCnvrtIntRef = ConvertSignalByDataset(sPreset.verticesPDF, mSigIntRef);


%         figure;
%         tiledlayout(2,3)
%         nexttile; imagesc(mSig); colorbar;title('Reference n')
%         nexttile; imagesc(mSigRecRep); colorbar; title('Rep rec')
%         nexttile; imagesc(mSigRecPhi); colorbar;title('ours rec')
%
%         nexttile; imagesc(mSigRef); colorbar;title('Reference N')
%         nexttile; imagesc(mSigRep); colorbar;title('Rep interp')
%         nexttile; imagesc(mSigInt); colorbar;title('ours interp')
%
%         figure;
%         tiledlayout(2,3)
%         nexttile; plot(mSig(:,1)); title('Reference n')
%         nexttile; plot(mSigRecRep(:,1)); title('Rep rec')
%         nexttile; plot(mSigRecPhi(:,1)); title('ours rec')
%
%         nexttile; plot(mSigRef(:,1)); title('Reference N')
%         nexttile; plot(mSigRep(:,1)); title('Rep interp')
%         nexttile; plot(mSigInt(:,1)); title('ours interp')


% ------------------------------------------------------------------------------------------
% Plots
% ------------------------------------------------------------------------------------------
sigIndToPlot     = min(7, size(mSigCnvrt,2));
vSigCoeffsV      = mSigCoeffsV(:,sigIndToPlot);
vSigRefCoeffsPhi = mSigRefCoeffsPhi(:,sigIndToPlot);
vSigCoeffsPhi    = mSigCoeffsPhi(:,sigIndToPlot);
vAlpha           = mAlpha(:,sigIndToPlot);

vSig       = mSigCnvrt(:,sigIndToPlot);
vSigRef    = mSigCnvrtRef(:,sigIndToPlot);
vSigRecV   = mSigCnvrtRecV(:,sigIndToPlot);
vSigRecPhi = mSigCnvrtRecPhi(:,sigIndToPlot);
vSigRecRep = mSigCnvrtRecRep(:,sigIndToPlot);
vSigInt    = mSigCnvrtInt(:,sigIndToPlot);
vSigIntRef = mSigCnvrtIntRef(:,sigIndToPlot);
vSigNys    = mSigCnvrtNys(:,sigIndToPlot);
vSigRep    = mSigCnvrtRep(:,sigIndToPlot);

if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
    figure; subplot(121); plot(mSigCoeffsPhi); subplot(122); plot(mAlpha)

    CRep = diag(lambdaPhi)*Phi.'*mAlpha;
    PlotCoeffsMatrix(mSigCoeffsPhi, '${\bf C}$', CRep, '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
    PlotCoefficients([], mSigCoeffsPhi(:,sigIndToPlot), lambdaPhi)
    PlotCoefficients(sPlotParams, vSigCoeffsPhi, lambdaPhi)
    PlotCoefficients(sPlotParams, vSigRefCoeffsPhi, lambdaPhi)
end


if sPlotParams.b_globalPlotEnable && sPreset.dim <=3
    if sPlotParams.b_plotExtraGraphSigAnalysis
        PlotExtraGraphSignalAnalysis(vSig, vSigRecPhi, vSigRecV, vSigRecRep, vSigRef, vSigInt, vSigNys, vSigRep,...
            vSigCoeffsPhi, vSigCoeffsV, vSigRefCoeffsPhi)
    end

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
    PlotTwoMoonsEigsRLS(sPlotParams, sDataset, sKernelParams, mSigRecPhi, mSigCoeffsPhi, sPreset.gamma1, sPreset.gamma2, sPreset.b_normalizePhi);
    PlotTwoMoonsLapRLS(sPlotParams, sDataset, sPreset.sDistanceParams, sPreset.omega, mAlpha, sPreset.gamma1Rep, sPreset.gamma2Rep);

elseif ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})
    b_transpose = strcmp(sPreset.verticesPDF, 'MNIST');
    nDigitsToPlot = 50;    
    if sPreset.n + nDigitsToPlot < sPreset.N
        vSamples = round(linspace(sPreset.n+1,sPreset.N,nDigitsToPlot));
        figTitle = ['Prediction on unseen $', num2str(length(vSamples)), '$ images'];
        figName = 'interpolation_signals';
        PlotDigits(sPlotParams, xInt(vSamples,:), mSigCnvrtInt(vSamples)-b_transpose, b_transpose, figTitle, figName);
    else
        vSamples = round(linspace(1,sPreset.n,nDigitsToPlot));
        figTitle = ['Prediction on seen $', num2str(length(vSamples)), '$ images'];
        PlotDigits([], xInt(vSamples,:), mSigCnvrtInt(vSamples)-b_transpose, b_transpose, figTitle);
    end

    vBadPredictInd = find(mSigCnvrtInt ~= mSigCnvrtRef);
    nBadPredict = min(50,length(vBadPredictInd));
    figTitle = ['Wrong prediction on given+unseen $n = ', num2str(nBadPredict), '$ points'];
    PlotDigits([], xInt(vBadPredictInd(1:nBadPredict),:), mSigCnvrtInt(vBadPredictInd(1:nBadPredict))-b_transpose, b_transpose, figTitle);

    nGmmPoints = 50;
    [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
    PhiGmm = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xGmm, sPreset.b_normalizePhi);
    mSigGmm = PhiGmm*mSigCoeffsPhi;
    mSigCnvrtGmm = ConvertSignalByDataset(sPreset.verticesPDF, mSigGmm);
    figTitle = [num2str(nGmmPoints), ' generated points'];
    figName = 'generated_signals';
    PlotDigits(sPlotParams, xGmm, mSigCnvrtGmm-b_transpose, b_transpose, figTitle, figName);

end
end