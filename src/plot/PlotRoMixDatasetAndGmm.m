function PlotRoMixDatasetAndGmm(sPlotParams, sPreset, sDataset, sDistParams)
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotDataVsGmm 
    nGmmPoints = 50*ismember(sPreset.verticesPDF, {'USPS', 'MNIST'}) + sPreset.n*~ismember(sPreset.verticesPDF, {'USPS', 'MNIST'});
    pltTitle = ['Dataset with n = ', num2str(sPreset.n), ' points'];
    plt2Title = ['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(sPreset.gmmNumComponents)];
    windowStyle = 'normal';
    PlotDataset(sPlotParams, sPreset, sDataset.sData.x, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle);
    if ismember(sPreset.verticesPDF, {'SwissRoll'}) && sPreset.sDatasetParams.b_randn
        b_transform = true;
        PlotDataset(sPlotParams, sPreset, sDataset.sData.S, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, b_transform);
    end
end
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotClustersAnalysis && sPreset.dim > 1
    %PlotGaussianEllipses(sPlotParams, sDistParams);
    b_plotCovMean = true; 
    PlotGmmResultWithDataset(sPlotParams, sDataset.sData.x, sDistParams, b_plotCovMean);
    b_plotCovMean = false; 
    PlotGmmResultWithDataset(sPlotParams, sDataset.sData.x, sDistParams, b_plotCovMean);
    PlotCovEigs(sPlotParams, sDistParams);
    if sDistParams.GMModel.NumComponents < 100
        PlotClustersMeans(sPreset, sPlotParams, sDistParams);
    end
end
end