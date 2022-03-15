function SaveFigure(sPlotParams, fig, figName, cFormatTypes)
if isfield(sPlotParams, 'outputFolder')
    if isfield(sPlotParams, 'actualDataDist') && isfield(sPlotParams, 'dim')
        simPrefix = strcat(sPlotParams.actualDataDist, '_', num2str(sPlotParams.dim), 'd_');
    else
        simPrefix = '';
    end
    for iFormat = 1:numel(cFormatTypes)
        format = cFormatTypes{iFormat};
        outputFolder = fullfile(sPlotParams.outputFolder, format);
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder)
        end
        saveas(fig, strcat(outputFolder, filesep, simPrefix, figName), format);
    end
end
end