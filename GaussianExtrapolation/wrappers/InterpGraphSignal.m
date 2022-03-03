function [mSigCnvrtRecPhi, mSigCnvrtRecV, mSigCnvrtRecRep, mSigCnvrt, mSigCnvrtInt, mSigCnvrtNys, mSigCnvrtRep, mSigCnvrtRef] = ...
    InterpGraphSignal(sPlotParams, sPreset, sDataset, sKernelParams, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, W, WRef, D, DRef, Ln, LnRef)
dim                = sPreset.dim;
n                  = sPreset.n;
N                  = sPreset.N;
MTilde             = sPreset.MTilde;
omega              = sPreset.omega; % for W
gamma1             = sPreset.gamma1;
gamma2             = sPreset.gamma2;
gamma1Rep          = sPreset.gamma1Rep;
gamma2Rep          = sPreset.gamma2Rep;
sDistanceParams    = sPreset.sDistanceParams;
b_normalizePhi     = sPreset.b_normalizePhi;
interpRatio        = N/n;

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
    mSigCoeffsPhi = EigsRLS(Phi, gamma1, gamma2, invLambda, Ln, mSigMasked, sPreset.b_maskDataFitTerm);
    b_normalizeAlpha = false;
    mAlpha = LapRLS(W, mSigMasked, Ln, gamma1Rep, gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
else
    mSigCoeffsPhi = EigsRLS(Phi, gamma1, gamma2, invLambda, Ln, mSig, sPreset.b_maskDataFitTerm);
    b_normalizeAlpha = false;
    mAlpha = LapRLS(W, mSig, Ln, gamma1Rep, gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
end

% Just for reference
mSigRefCoeffsPhi = EigsRLS(PhiInt, gamma1, gamma2, invLambda, LnRef, mSigRef, sPreset.b_maskDataFitTerm);
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


if sPlotParams.b_globalPlotEnable && dim <=3
    if sPlotParams.b_plotExtraGraphSigAnalysis
        PlotExtraGraphSignalAnalysis(vSig, vSigRecPhi, vSigRecV, vSigRecRep, vSigRef, vSigInt, vSigNys, vSigRep,...
            vSigCoeffsPhi, vSigCoeffsV, vSigRefCoeffsPhi)
    end

    if dim == 1
        colorOrder = get(gca, 'ColorOrder');
        PlotGraphSignals(sPlotParams, ['Graph signal reconstruction'], 'Reconstruction', ...
            {xTrain, xTrain}, ...
            {vSig, vSigRecPhi}, ...
            {'$s$', '$s_{\Phi}^{{\bf rec}}$'}, ...
            {n, n}, {'o', '.'}, mat2cell(colorOrder(1:2,:),[1 1], 3));
        PlotGraphSignals(sPlotParams, ['Graph signal interpolation'], 'Interpolation', ...
            {xInt, xInt}, ...
            {vSigRef, vSigInt}, ...
            {'$s^{{\bf ref}}$', '$s^{{\bf int}}$'}, ...
            {n, n}, {'o', '.'}, mat2cell(colorOrder(3:4,:),[1 1], 3));
    else
        PlotGraphSignals(sPlotParams, ['Graph signal interpolation'], 'Interpolation', ...
            {xTrain, xInt, xTrain, xInt}, ...
            {vSig, vSigRef, vSigRecPhi, vSigInt}, ...
            {'$s$', '$s^{{\bf ref}}$', '$s_{\Phi}^{{\bf rec}}$', '$s^{{\bf int}}$'}, ...
            {n, n, n, n});
    end
%     %     PlotGraphSignals(sPlotParams, ['Graph signals on given $n=' num2str(n) '$ nodes (Train set)'], 'TrainSet', ...
%     %         {xTrain, xTrain, xTrain, xTrain}, ...
%     %         {vSig, vSigRecPhi, vSigRecV, vSigRecRep}, ...
%     %         {'$s$', '$s_{\Phi}^{{\bf rec}}$', '$s_{V}^{{\bf rec}}$', '$s_{K}^{{\bf rec}}$'}, ...
%     %         {n, n, n, n});
%     %     cmap = PlotGraphSignals(sPlotParams, ['Graph signals on all $N=' num2str(N) '$ nodes (Train \& Test sets)'], 'TrainAndTestSet', ...
%     %         {xInt, xInt, xInt, xInt}, ...
%     %         {vSigRef, vSigInt, vSigNys, vSigRep}, ...
%     %         {'$s^{{\bf ref}}$', '$s^{{\bf int}}$', '$s^{{\bf nys}}$', '$s^{{\bf rep}}$'}, ...
%     %         {n, n, n, n});
    if sPlotParams.b_plotGmmSignal
        % GMM
        nGmmPoints = 1000;
        [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
        [PhiGmm, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xGmm, b_normalizePhi);
        mSigGmm = PhiGmm*mSigCoeffsPhi;
        vSigGmm = mSigGmm(:,sigIndToPlot);
        if dim == 2
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
    PlotTwoMoonsEigsRLS(sPlotParams, sDataset, sKernelParams, mSigRecPhi, mSigCoeffsPhi, gamma1,gamma2, b_normalizePhi);
    PlotTwoMoonsLapRLS(sPlotParams, sDataset, sDistanceParams, omega, mAlpha, gamma1Rep, gamma2Rep);

elseif ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})
    b_transpose = strcmp(sPreset.verticesPDF, 'MNIST');
    nDigitsToPlot = 50;    
    if sPreset.n + nDigitsToPlot < sPreset.N
        vSamples = round(linspace(sPreset.n+1,sPreset.N,nDigitsToPlot));
        figTitle = ['Prediction on unseen $', num2str(length(vSamples)), '$ images'];
        figName = 'interpolation_signals';
        PlotDigits(sPlotParams, xInt(vSamples,:), mSigCnvrtInt(vSamples)-b_transpose, b_transpose, figTitle, figName)
    else
        vSamples = round(linspace(1,sPreset.n,nDigitsToPlot));
        figTitle = ['Prediction on seen $', num2str(length(vSamples)), '$ images'];
        PlotDigits([], xInt(vSamples,:), mSigCnvrtInt(vSamples)-b_transpose, b_transpose, figTitle)
    end

    vBadPredictInd = find(mSigCnvrtInt ~= mSigCnvrtRef);
    nBadPredict = min(50,length(vBadPredictInd));
    figTitle = ['Wrong prediction on given+unseen $n = ', num2str(nBadPredict), '$ points'];
    PlotDigits([], xInt(vBadPredictInd(1:nBadPredict),:), mSigCnvrtInt(vBadPredictInd(1:nBadPredict))-b_transpose, b_transpose, figTitle)

    nGmmPoints = 50;
    [xGmm,compIdx] = random(sKernelParams.sDistParams.GMModel, nGmmPoints);
    PhiGmm = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xGmm, b_normalizePhi);
    mSigGmm = PhiGmm*mSigCoeffsPhi;
    mSigCnvrtGmm = ConvertSignalByDataset(sPreset.verticesPDF, mSigGmm);
    figTitle = [num2str(nGmmPoints), ' generated points'];
    figName = 'generated_signals';
    PlotDigits(sPlotParams, xGmm, mSigCnvrtGmm-b_transpose, b_transpose, figTitle, figName)

end
end