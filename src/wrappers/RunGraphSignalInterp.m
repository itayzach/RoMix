function sResults = RunGraphSignalInterp(sPreset, sPlotParams, b_interpEigenvecs)

vNumLabeled = sPreset.nLabeled;
for iRun = 1:numel(vNumLabeled)
    sPreset = UpdatePreset(sPreset, b_interpEigenvecs, vNumLabeled(iRun));
    nMethods = numel(sPreset.cMethods);
    tSigCnvrtRec    = zeros(sPreset.n, sPreset.nSignals, sPreset.R, nMethods);
    tSigCnvrtInt    = zeros(sPreset.N, sPreset.nSignals, sPreset.R, nMethods);
    tSigCnvrtRecRef = zeros(sPreset.n, sPreset.nSignals, sPreset.R);
    tSigCnvrtIntRef = zeros(sPreset.N, sPreset.nSignals, sPreset.R);
    [tTrainTime, tIntTime] = deal(zeros(sPreset.R, nMethods));
    for r = 1:sPreset.R
        % Plot only on first iteration
        sPlotParams.b_globalPlotEnable = (iRun == 1) && (r == 1) && sPlotParams.b_globalPlotEnable;
        
        % Generate dataset
        fprintf('Run %d/%d. nLabeled = %d. Iteration r = %d of R = %d\n',iRun,numel(vNumLabeled),vNumLabeled(iRun),r,sPreset.R)
        t = tic;
        sDataset = GenerateDataset(sPlotParams, sPreset, b_interpEigenvecs);
    
        % Run
        for methodInd = 1:numel(sPreset.cMethods)
            funcHandle = GetMethodHandleFromStr(sPreset.cMethods{methodInd});
            [tSigCnvrtRec(:,:,r,methodInd), tSigCnvrtInt(:,:,r,methodInd), tTrainTime(r,methodInd), tIntTime(r,methodInd)] = funcHandle(sPlotParams, sPreset, sDataset);
        end
        
        % Convert expected signals
        tSigCnvrtRecRef(:,:,r) = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);
        tSigCnvrtIntRef(:,:,r) = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.yt);
        fprintf('Iteration r = %d of R = %d finished (took %.2f min)\n',r,sPreset.R,toc(t)/60)
    end
    [mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vTrainTimeStd, vIntTime, vIntTimeStd] = ...
        InterpGraphSignalsMetrics(sPlotParams, sPreset, b_interpEigenvecs, tSigCnvrtRec, tSigCnvrtRecRef, tSigCnvrtInt, tSigCnvrtIntRef, tTrainTime, tIntTime);
    sResults(iRun) = SaveResults(sPreset, sPlotParams, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vTrainTimeStd, vIntTime, vIntTimeStd);
end
if numel(vNumLabeled) > 1
    PlotMetricsVsNumLabeled(sPreset, sPlotParams, sResults)
end
end