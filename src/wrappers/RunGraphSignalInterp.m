function sResults = RunGraphSignalInterp(sPreset, sPlotParams, b_interpEigenvecs)

sPreset = UpdatePreset(sPreset, b_interpEigenvecs);
nMethods = numel(sPreset.cMethods);

tSigCnvrtRec    = zeros(sPreset.n, sPreset.nSignals, sPreset.R, nMethods);
tSigCnvrtInt    = zeros(sPreset.N, sPreset.nSignals, sPreset.R, nMethods);
tSigCnvrtRecRef = zeros(sPreset.n, sPreset.nSignals, sPreset.R);
tSigCnvrtIntRef = zeros(sPreset.N, sPreset.nSignals, sPreset.R);
[tTrainTime, tIntTime] = deal(zeros(sPreset.R, nMethods));
for r = 1:sPreset.R
    % Plot only on first iteration
    sPlotParams.b_globalPlotEnable = (r == 1) && sPlotParams.b_globalPlotEnable;
    
    % Generate dataset
    fprintf('Iteration r = %d of R = %d\n',r,sPreset.R)
    t = tic;
    sDataset = GenerateDataset(sPlotParams, sPreset, b_interpEigenvecs);
    tSigCnvrtRecRef(:,:,r) = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
    tSigCnvrtIntRef(:,:,r) = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);

    % Run
    [tSigCnvrtRec(:,:,r,1), tSigCnvrtInt(:,:,r,1), tTrainTime(r,1), tIntTime(r,1)] = InterpGraphSignalRoMix(sPlotParams, sPreset, sDataset);
    if sPreset.b_compareMethods
        [tSigCnvrtRec(:,:,r,2), tSigCnvrtInt(:,:,r,2), tTrainTime(r,2), tIntTime(r,2)] = InterpGraphSignalNystrom(sPlotParams, sPreset, sDataset);
        [tSigCnvrtRec(:,:,r,3), tSigCnvrtInt(:,:,r,3), tTrainTime(r,3), tIntTime(r,3)] = InterpGraphSignalRepThm(sPlotParams, sPreset, sDataset);
        [tSigCnvrtRec(:,:,r,4), tSigCnvrtInt(:,:,r,4), tTrainTime(r,4), tIntTime(r,4)] = InterpGraphSignalPesenson(sPlotParams, sPreset, sDataset);
        [tSigCnvrtRec(:,:,r,5), tSigCnvrtInt(:,:,r,5), tTrainTime(r,5), tIntTime(r,5)] = InterpGraphSignalKnn(sPlotParams, sPreset, sDataset);
    end
    fprintf('Iteration r = %d of R = %d finished (took %.2f min)\n',r,sPreset.R,toc(t)/60)
end

[mAccRec, mAccInt] = InterpGraphSignalsMetrics(sPlotParams, sPreset, b_interpEigenvecs, tSigCnvrtRec, tSigCnvrtRecRef, tSigCnvrtInt, tSigCnvrtIntRef, tTrainTime, tIntTime);
sResults = SaveResults(sPreset, sPlotParams, tSigCnvrtRec, tSigCnvrtInt, tSigCnvrtRecRef, tSigCnvrtIntRef, mAccRec, mAccInt, tTrainTime, tIntTime);

end