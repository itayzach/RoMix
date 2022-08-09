function InterpEigenvecsMetrics(sPlotParams, sPreset, ...
    tVIntToCompare, tVRefToCompare, tVRecPhiToCompare, tVToCompare, tVNysToCompare, tVRepToCompare)
[vAccInt, vAccStdInt] = CalcErrAndAcc(tVIntToCompare, tVRefToCompare, 'Analytic');
[vAccRecPhi, vAccStdPhi] = CalcErrAndAcc(tVRecPhiToCompare, tVToCompare, 'Analytic');
if sPreset.b_compareMethods
    [vAccNys, vAccStdNys] = CalcErrAndAcc(tVNysToCompare, tVRefToCompare, 'Nystrom');
    [vAccRep, vAccStdRep] = CalcErrAndAcc(tVRepToCompare, tVRefToCompare, 'Representer');
    PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], [vAccStdInt, vAccStdNys, vAccStdRep], ...
        {'Acc$(\tilde{\psi}^{{\bf RoMix}}_m, \tilde{\psi}_m)$', 'Acc$(\tilde{\psi}^{{\bf nys}}_m, \tilde{\psi}_m)$', ...
        'Acc$(\tilde{\psi}^{{\bf rep}}_m, \tilde{\psi}_m)$'}, [ sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)]);
else
    PlotAccuracy(sPlotParams, [vAccInt, vAccRecPhi], [vAccStdInt vAccStdPhi], ...
        {'Acc$(\tilde{\psi}^{{\bf RoMix}}_m, \tilde{\psi}_m)$', 'Acc$(\psi^{{\bf RoMix}}_m, \psi_m)$'}, [sPreset.matrixForEigs '_Acc_eigs_0_to_' num2str(sPreset.M-1)]);
end
end