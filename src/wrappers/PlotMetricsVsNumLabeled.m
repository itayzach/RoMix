function PlotMetricsVsNumLabeled(sPreset, sPlotParams, sResults)
b_zoomIn = sPreset.sDatasetParams.b_zoomInAcc;
xTickNames = sPreset.sDatasetParams.xTickNames;

mAccRec = reshape([sResults(:).mAccRec], numel(sPreset.cMethods),[])';
mAccStdRec = reshape([sResults(:).mAccStdRec], numel(sPreset.cMethods),[])';
PlotMetric(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ...
    'nLabelInterpAcc', [], xTickNames, 'Interpolation accuracy', [], [], b_zoomIn);

mAccInt = reshape([sResults(:).mAccInt], numel(sPreset.cMethods),[])';
mAccStdInt = reshape([sResults(:).mAccStdInt], numel(sPreset.cMethods),[])';
PlotMetric(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
    'nLabelExtrapAcc', [], xTickNames, 'Extrapolation accuracy', [], [], b_zoomIn);

mInterpTime = reshape([sResults(:).vTrainTime], numel(sPreset.cMethods),[])';
mInterpTimeStd = reshape([sResults(:).vTrainTimeStd], numel(sPreset.cMethods),[])';
PlotMetric(sPlotParams, mInterpTime, mInterpTimeStd, sPreset.cMethods, ...
    'nLabelInterpTime', [], xTickNames, 'Interpolation time', 'Time [s]');

mExtrapTime = reshape([sResults(:).vIntTime], numel(sPreset.cMethods),[])';
mExtrapTimeStd = reshape([sResults(:).vIntTimeStd], numel(sPreset.cMethods),[])';
PlotMetric(sPlotParams, mExtrapTime, mExtrapTimeStd, sPreset.cMethods, ...
    'nLabelExtrapTime', [], xTickNames, 'Extrapolation time', 'Time [s]');
end