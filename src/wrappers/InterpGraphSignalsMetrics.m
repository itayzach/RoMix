function [mAccRec, mAccInt, vTrainTime, vIntTime] = ...
    InterpGraphSignalsMetrics(sPlotParams, sPreset, b_interpEigenvecs, tSigCnvrtRec, tSigCnvrtRecRef, tSigCnvrtInt, tSigCnvrtIntRef, tTrainTime, tIntTime)
[mAccRec(:,1), mAccStdRec(:,1)] = CalcErrAndAcc(squeeze(tSigCnvrtRec(:,:,:,1)), tSigCnvrtRecRef);
[mAccInt(:,1), mAccStdInt(:,1)] = CalcErrAndAcc(squeeze(tSigCnvrtInt(:,:,:,1)), tSigCnvrtIntRef);
vTrainTime = mean(tTrainTime,1);
vIntTime = mean(tIntTime,1);
if sPreset.b_compareMethods
    for methodInd = 2:size(tSigCnvrtRec,4)
        [mAccRec(:,methodInd), mAccStdRec(:,methodInd)] = CalcErrAndAcc(tSigCnvrtRec(:,:,:,methodInd), tSigCnvrtRecRef);
    end
    for methodInd = 2:size(tSigCnvrtInt,4)
        [mAccInt(:,methodInd), mAccStdInt(:,methodInd)] = CalcErrAndAcc(tSigCnvrtInt(:,:,:,methodInd), tSigCnvrtIntRef);
    end
    if sPreset.nSignals == 1
        PrintAccuracyLatex(sPreset, mAccRec, mAccStdRec, mAccInt, mAccStdInt, sPreset.cMethods);
        PrintAccuracy(sPreset, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vIntTime, sPreset.cMethods);
    elseif b_interpEigenvecs
        PlotAccuracy(sPlotParams, mAccInt, mAccStdInt, ...
            sPreset.cMethods, [ sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)]);
    end
    if isfield(sPreset.sDatasetParams, 'monthNames')
        PlotAccuracy(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
            'InterpAcc', [], sPreset.sDatasetParams.monthNames, 'Interpolation accuracy');
        %{'Acc$(\tilde{s}^{{\bf RoMix}}_m, \tilde{s}_m)$', 'Acc$(\tilde{s}^{{\bf nys}}_m, \tilde{s}_m)$', 'Acc$(\tilde{s}^{{\bf rep}}_m, \tilde{s}_m)$', 'Acc$(\tilde{s}^{{\bf PW}}_m, \tilde{s}_m)$', 'Acc$(\tilde{s}^{{\bf kNN}}_m, \tilde{s}_m)$'}
        PlotAccuracy(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ... 
            'ProjAcc', [], sPreset.sDatasetParams.monthNames, 'Projection accuracy');
            %{'Acc$(s^{{\bf RoMix}}_m, s_m)$', 'Acc$(s^{{\bf nys}}_m, s_m)$', 'Acc$(s^{{\bf rep}}_m, s_m)$', 'Acc$(s^{{\bf PW}}_m, s_m)$', 'Acc$(s^{{\bf kNN}}_m, s_m)$'}, ...
    end
else
    if sPreset.nSignals == 1
        PrintAccuracyLatex(sPreset, mAccRec(:,1), mAccStdRec(:,1), mAccInt(:,1), mAccStdInt(:,1), {'RoMix'});
        PrintAccuracy(sPreset, mAccRec(:,1), mAccStdRec(:,1), mAccInt(:,1), mAccStdInt(:,1), vTrainTime, vIntTime, {'RoMix'});
    elseif b_interpEigenvecs
        PlotAccuracy(sPlotParams, [mAccInt(:,1), mAccRec(:,1)], [mAccStdInt(:,1), mAccStdRec(:,1)], ...
            {'Acc$(\tilde{\psi}^{{\bf RoMix}}_m, \tilde{\psi}_m)$', 'Acc$(\psi^{{\bf RoMix}}_m, \psi_m)$'}, ...
            [sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)]);
    end
    if isfield(sPreset.sDatasetParams, 'monthNames')
        PlotAccuracy(sPlotParams, [mAccInt(:,1), mAccRec(:,1)], [mAccStdInt(:,1) mAccStdRec(:,1)], ...
            {'Acc$(\tilde{s}^{{\bf RoMix}}_m, \tilde{s}_m)$', 'Acc$(s^{{\bf RoMix}}_m, s_m)$'}, ...
            'ProjInterpAcc', [], sPreset.sDatasetParams.monthNames, 'Projection \& interpolation accuracy');
    end
end
end