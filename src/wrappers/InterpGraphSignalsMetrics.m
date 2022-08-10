function [mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vTrainTimeStd, vIntTime, vIntTimeStd] = ...
    InterpGraphSignalsMetrics(sPlotParams, sPreset, b_interpEigenvecs, tSigCnvrtRec, tSigCnvrtRecRef, tSigCnvrtInt, tSigCnvrtIntRef, tTrainTime, tIntTime)
[mAccRec(:,1), mAccStdRec(:,1)] = CalcErrAndAcc(tSigCnvrtRec(:,:,:,1), tSigCnvrtRecRef);
[mAccInt(:,1), mAccStdInt(:,1)] = CalcErrAndAcc(tSigCnvrtInt(:,:,:,1), tSigCnvrtIntRef);
vTrainTime = mean(tTrainTime,1);
vIntTime = mean(tIntTime,1);
vTrainTimeStd = std(tTrainTime,1);
vIntTimeStd = std(tIntTime,1);
if sPreset.b_compareMethods
    for methodInd = 2:size(tSigCnvrtRec,4)
        [mAccRec(:,methodInd), mAccStdRec(:,methodInd)] = CalcErrAndAcc(tSigCnvrtRec(:,:,:,methodInd), tSigCnvrtRecRef);
    end
    for methodInd = 2:size(tSigCnvrtInt,4)
        [mAccInt(:,methodInd), mAccStdInt(:,methodInd)] = CalcErrAndAcc(tSigCnvrtInt(:,:,:,methodInd), tSigCnvrtIntRef);
    end
    if sPreset.nSignals == 1
        PrintAccuracyLatex(sPreset, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, sPreset.cMethods);
        PrintAccuracy(sPreset, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vIntTime, sPreset.cMethods);
    elseif b_interpEigenvecs
        PlotAccuracy(sPlotParams, mAccInt, mAccStdInt, ...
            sPreset.cMethods, [ sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)], [95 100]);
    end
    if isfield(sPreset.sDatasetParams, 'xTickNames') && sPreset.nSignals > 1
        PlotAccuracy(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
            'ExtrapAcc', [], sPreset.sDatasetParams.xTickNames, 'Extrapolation accuracy');
        PlotAccuracy(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ... 
            'InterpAcc', [], sPreset.sDatasetParams.xTickNames, 'Interpolation accuracy');
    end
else
    if sPreset.nSignals == 1
        PrintAccuracyLatex(sPreset, mAccRec(:,1), mAccStdRec(:,1), mAccInt(:,1), mAccStdInt(:,1), vTrainTime(1), {'RoMix'});
        PrintAccuracy(sPreset, mAccRec(:,1), mAccStdRec(:,1), mAccInt(:,1), mAccStdInt(:,1), vTrainTime(1), vIntTime(1), {'RoMix'});
    elseif b_interpEigenvecs
        PlotAccuracy(sPlotParams, [mAccInt(:,1), mAccRec(:,1)], [mAccStdInt(:,1), mAccStdRec(:,1)], ...
            {'Acc$(\tilde{\psi}^{{\bf RoMix}}_m, \tilde{\psi}_m)$', 'Acc$(\psi^{{\bf RoMix}}_m, \psi_m)$'}, ...
            [sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)], [95 100]);
    end
    if isfield(sPreset.sDatasetParams, 'xTickNames')
        PlotAccuracy(sPlotParams, [mAccInt(:,1), mAccRec(:,1)], [mAccStdInt(:,1) mAccStdRec(:,1)], ...
            {'Acc$(\tilde{s}^{{\bf RoMix}}_m, \tilde{s}_m)$', 'Acc$(s^{{\bf RoMix}}_m, s_m)$'}, ...
            'InterpExtrapAcc', [], sPreset.sDatasetParams.xTickNames, 'Interpolation \& extrapolation accuracy');
    end
end
end