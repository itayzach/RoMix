function PlotClustersMeans(sPreset, sDistParams)
if ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})
    b_transpose = strcmp(sPreset.verticesPDF, 'MNIST');
    vSamples = 1:sDistParams.estNumComponents;
    xMeans = cell2mat(sDistParams.mu');
    PlotDigits([], xMeans(vSamples,:), vSamples, b_transpose, 'Means');
end
end