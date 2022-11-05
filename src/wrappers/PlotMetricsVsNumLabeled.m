function PlotMetricsVsNumLabeled(sPreset, sPlotParams, sResults)
b_zoomIn = sPreset.sDatasetParams.b_zoomInAcc;
b_compareDb = ismember(sPreset.verticesPDF, {'BulgariBeacons'});
xTickNames = sPreset.sDatasetParams.xTickNames;
if ismember(sPreset.verticesPDF, {'MNIST'})
    xlab = '$\ell / n [\%]$';
else
    xlab = '$\ell$';
end

mAccRec = reshape([sResults(:).mAccRec], numel(sPreset.cMethods),[])';
mAccStdRec = reshape([sResults(:).mAccStdRec], numel(sPreset.cMethods),[])';
mAccInt = reshape([sResults(:).mAccInt], numel(sPreset.cMethods),[])';
mAccStdInt = reshape([sResults(:).mAccStdInt], numel(sPreset.cMethods),[])';
if b_compareDb
    PlotMetric(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ...
        'nLabelInterpErr', [], xTickNames, [], 'Error [dB/node]', xlab, b_zoomIn);
    PlotMetric(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
        'nLabelExtrapErr', [], xTickNames, [], 'Error [dB/node]', xlab, b_zoomIn);
else
    if ismember(sPreset.verticesPDF, {'MNIST'})
        xylim = [0 100];
    else
        xylim = [];
    end
    PlotMetric(sPlotParams, mAccRec, mAccStdRec, sPreset.cMethods, ...
        'nLabelInterpAcc', xylim, xTickNames, [], 'Accuracy $[\%]$', xlab, b_zoomIn);
    PlotMetric(sPlotParams, mAccInt, mAccStdInt, sPreset.cMethods, ...
        'nLabelExtrapAcc', xylim, xTickNames, [], 'Accuracy $[\%]$', xlab, b_zoomIn);
end

b_zoomIn = false;
b_ylog = ismember(sPreset.verticesPDF, {'BulgariBeacons'});
mInterpTime = reshape([sResults(:).vTrainTime], numel(sPreset.cMethods),[])';
mInterpTimeStd = reshape([sResults(:).vTrainTimeStd], numel(sPreset.cMethods),[])';
PlotMetric(sPlotParams, mInterpTime, mInterpTimeStd, sPreset.cMethods, ...
    'nLabelInterpTime', [], xTickNames, [], 'Time [s]', xlab, b_zoomIn, b_ylog);

%b_zoomIn = false;
%b_ylog = ismember(sPreset.verticesPDF, {'BulgariBeacons'});
%mExtrapTime = reshape([sResults(:).vIntTime], numel(sPreset.cMethods),[])';
%mExtrapTimeStd = reshape([sResults(:).vIntTimeStd], numel(sPreset.cMethods),[])';
%PlotMetric(sPlotParams, mExtrapTime, mExtrapTimeStd, sPreset.cMethods, ...
%   'nLabelExtrapTime', [], xTickNames, 'Extrapolation time', 'Time [s]', xlab, b_zoomIn, b_ylog);
