function InterpGraphSignalsMetrics(sPlotParams, sPreset, ...
    tSigCnvrtRecPhi, tSigCnvrt, tSigCnvrtInt, tSigCnvrtRef, tSigCnvrtRecRep, tSigCnvrtRep, tSigCnvrtRecV, tSigCnvrtNys)
[vAccRecPhi, vAccStdPhi] = CalcErrAndAcc(tSigCnvrtRecPhi, tSigCnvrt, 'RoMix (train)');
[vAccInt, vAccStdInt]    = CalcErrAndAcc(tSigCnvrtInt, tSigCnvrtRef, 'RoMix (test)');
if sPreset.b_compareMethods
    [vAccRecRep, vAccStdRecRep] = CalcErrAndAcc(tSigCnvrtRecRep, tSigCnvrt, 'Representer (train)');
    [vAccRep, vAccStdRep]       = CalcErrAndAcc(tSigCnvrtRep, tSigCnvrtRef, 'Representer (test)');
    [vAccRecV, vAccStdRecV]     = CalcErrAndAcc(tSigCnvrtRecV, tSigCnvrt, 'Nystrom (train)');
    [vAccNys, vAccStdNys]       = CalcErrAndAcc(tSigCnvrtNys, tSigCnvrtRef, 'Nystrom (test)');
    if isscalar(vAccRecPhi)
        PrintAccuracyLatex(sPreset, [vAccRecPhi, vAccRecRep, vAccRecV], [vAccStdPhi, vAccStdRecRep, vAccStdRecV], ...
            [vAccInt, vAccRep, vAccNys], [vAccStdInt, vAccStdRep, vAccStdNys], {'RoMix', 'Rep. Thm', 'Nystrom'});
    end
    if isfield(sPreset.sDatasetParams, 'monthNames')
        PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], [vAccStdInt, vAccStdNys, vAccStdRep], ...
            {'Acc$(\tilde{s}^{{\bf RoMix}}_m, \tilde{s}_m)$', 'Acc$(\tilde{s}^{{\bf nys}}_m, \tilde{s}_m)$', ...
            'Acc$(\tilde{s}^{{\bf rep}}_m, \tilde{s}_m)$'}, 'InterpAcc', [], sPreset.sDatasetParams.monthNames, 'Interpolation accuracy');
        PlotAccuracy(sPlotParams, [vAccRecPhi, vAccRecV, vAccRecRep], [vAccStdPhi, vAccStdRecV, vAccStdRecRep], ...
            {'Acc$(s^{{\bf RoMix}}_m, s_m)$', 'Acc$(s^{{\bf nys}}_m, s_m)$', ...
            'Acc$(s^{{\bf rep}}_m, s_m)$'}, 'ProjAcc', [], sPreset.sDatasetParams.monthNames, 'Projection accuracy');
    end
else
    if isscalar(vAccRecPhi)
        PrintAccuracyLatex(sPreset, vAccRecPhi, vAccStdPhi, vAccInt, vAccStdInt, {'RoMix'});
    end
    if isfield(sPreset.sDatasetParams, 'monthNames')
        PlotAccuracy(sPlotParams, [vAccInt, vAccRecPhi], [vAccStdInt vAccStdPhi], ...
            {'Acc$(\tilde{s}^{{\bf RoMix}}_m, \tilde{s}_m)$', 'Acc$(s^{{\bf RoMix}}_m, s_m)$'}, ...
            'ProjInterpAcc', [], sPreset.sDatasetParams.monthNames, 'Projection \& interpolation accuracy');
    end
end
end