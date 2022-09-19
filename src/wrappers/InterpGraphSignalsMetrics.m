function [mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vTrainTimeStd, vIntTime, vIntTimeStd] = ...
    InterpGraphSignalsMetrics(sPlotParams, sPreset, b_interpEigenvecs, tSigCnvrtRec, tSigCnvrtRecRef, tSigCnvrtInt, tSigCnvrtIntRef, tSigCnvrtLabRef, tTrainTime, tIntTime)

assert(isequal(tSigCnvrtLabRef(sPreset.nLabeled+1:end,:,:), zeros(sPreset.n-sPreset.nLabeled, sPreset.nSignals, sPreset.R)));
assert(isequal(tSigCnvrtLabRef(1:sPreset.nLabeled,:,:), tSigCnvrtRecRef(1:sPreset.nLabeled,:,:)));
assert(isequal(tSigCnvrtLabRef(1:sPreset.nLabeled,:,:), tSigCnvrtIntRef(1:sPreset.nLabeled,:,:)));
assert(isequal(tSigCnvrtRecRef, tSigCnvrtIntRef(1:sPreset.n,:,:)));

vRecInd = 1:sPreset.n; %sPreset.nLabeled+1:sPreset.n;
vIntInd = 1:sPreset.N; %sPreset.n+1:sPreset.N;
b_compareDb = ismember(sPreset.verticesPDF, {'BulgariBeacons'});

[mAccRec(:,1), mAccStdRec(:,1)] = CalcErrAndAcc(tSigCnvrtRec(vRecInd,:,:,1), tSigCnvrtRecRef(vRecInd,:,:), b_compareDb);
[mAccInt(:,1), mAccStdInt(:,1)] = CalcErrAndAcc(tSigCnvrtInt(vIntInd,:,:,1), tSigCnvrtIntRef(vIntInd,:,:), b_compareDb);
vTrainTime = mean(tTrainTime,1);
vIntTime = mean(tIntTime,1);
vTrainTimeStd = std(tTrainTime,[],1);
vIntTimeStd = std(tIntTime,[],1);
if sPreset.b_compareMethods
    for methodInd = 2:size(tSigCnvrtRec,4)
        [mAccRec(:,methodInd), mAccStdRec(:,methodInd)] = CalcErrAndAcc(tSigCnvrtRec(vRecInd,:,:,methodInd), tSigCnvrtRecRef(vRecInd,:,:), b_compareDb);
    end
    for methodInd = 2:size(tSigCnvrtInt,4)
        [mAccInt(:,methodInd), mAccStdInt(:,methodInd)] = CalcErrAndAcc(tSigCnvrtInt(vIntInd,:,:,methodInd), tSigCnvrtIntRef(vIntInd,:,:), b_compareDb);
    end
    if sPreset.nSignals == 1
        PrintAccuracyLatex(sPreset, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, sPreset.cMethods);
        PrintAccuracy(sPreset, mAccRec, mAccStdRec, mAccInt, mAccStdInt, vTrainTime, vIntTime, sPreset.cMethods);
    elseif b_interpEigenvecs
        PlotMetric(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
            [ sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)], [97 100]);
    end
    if isfield(sPreset.sDatasetParams, 'xTickNames') && sPreset.nSignals > 1
        PlotMetric(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
            'ExtrapAcc', [], sPreset.sDatasetParams.xTickNames, 'Extrapolation accuracy');
        PlotMetric(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ... 
            'InterpAcc', [], sPreset.sDatasetParams.xTickNames, 'Interpolation accuracy');
    end
else
    if sPreset.nSignals == 1
        PrintAccuracyLatex(sPreset, mAccRec(:,1), mAccStdRec(:,1), mAccInt(:,1), mAccStdInt(:,1), vTrainTime(1), {'RoMix'});
        PrintAccuracy(sPreset, mAccRec(:,1), mAccStdRec(:,1), mAccInt(:,1), mAccStdInt(:,1), vTrainTime(1), vIntTime(1), {'RoMix'});
    elseif b_interpEigenvecs
        PlotMetric(sPlotParams, [mAccRec(:,1), mAccInt(:,1)], [mAccStdRec(:,1), mAccStdInt(:,1)], ...
            {'Interp.', 'Extrap.'}, ... {'Acc$(\psi^{{\bf RoMix}}_m, \psi_m)$', 'Acc$(\tilde{\psi}^{{\bf RoMix}}_m, \tilde{\psi}_m)$'}, ...
            [sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)], [97 100]);
    end
    if isfield(sPreset.sDatasetParams, 'xTickNames')
        PlotMetric(sPlotParams, [mAccRec(:,1), mAccInt(:,1)], [mAccStdRec(:,1), mAccStdInt(:,1)], ...
            {'Interp.', 'Extrap.'}, ...{'Acc$(s^{{\bf RoMix}}_m, s_m)$', 'Acc$(\tilde{s}^{{\bf RoMix}}_m, \tilde{s}_m)$'}, ...
            'InterpExtrapAcc', [], sPreset.sDatasetParams.xTickNames, 'Interpolation \& extrapolation accuracy');

    end
end
end