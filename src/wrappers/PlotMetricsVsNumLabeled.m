function PlotMetricsVsNumLabeled(sPreset, sPlotParams, sResults)
b_zoomIn = sPreset.sDatasetParams.b_zoomInAcc;
xTickNames = sPreset.sDatasetParams.xTickNames;

mAccRec = reshape([sResults(:).mAccRec], numel(sPreset.cMethods),[])';
mAccStdRec = reshape([sResults(:).mAccStdRec], numel(sPreset.cMethods),[])';
PlotAccuracy(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ...
    'nLabelInterpAcc', [], xTickNames, 'Interpolation accuracy', [], [], b_zoomIn);

mAccInt = reshape([sResults(:).mAccInt], numel(sPreset.cMethods),[])';
mAccStdInt = reshape([sResults(:).mAccStdInt], numel(sPreset.cMethods),[])';
PlotAccuracy(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
    'nLabelExtrapAcc', [], xTickNames, 'Extrapolation accuracy', [], [], b_zoomIn);

mTrainTime = reshape([sResults(:).vTrainTime], numel(sPreset.cMethods),[])';
mTrainTimeStd = reshape([sResults(:).vTrainTimeStd], numel(sPreset.cMethods),[])';
PlotAccuracy(sPlotParams, mTrainTime, mTrainTimeStd, sPreset.cMethods, ...
    'nLabelTrainTime', [], xTickNames, 'Interpolation time', 'Time [s]');
end