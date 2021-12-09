function PlotGraphSignalErrors(sPlotParams, xInd, sig1, cSignals, cDispName, pltTitle)

cSignals = reshape(cSignals,[],1);
cDispName = reshape(cDispName,[],1);
n = length(sig1);
cMarkerType = {'x', '+', '*'};
fig = figure('Name', 'Graph signals errors');
for sigInd = 1:numel(cSignals)
    sigDispName = cDispName{sigInd};
    sig = cSignals{sigInd};
    plot(xInd(1):xInd(2), abs(sig1-sig), cMarkerType{sigInd}, 'DisplayName', ['$' sigDispName '$']);
    hold on;
end
for sigInd = 1:numel(cSignals)
    sigDispName = cDispName{sigInd};
    sig = cSignals{sigInd};
    plot(xInd(1):xInd(2), mean(abs(sig1-sig))*ones(1,n), 'LineWidth', 3, 'DisplayName', ['${{\bf avg}}({' sigDispName '})$']);
    hold on;
end
xlim([xInd(1),xInd(2)]);
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title(pltTitle, 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
x0     = 10;
y0     = 50;
height = 350;
width  = 600;
set(gcf,'Position', [x0 y0 width height])

%% Save
if isfield(sPlotParams, 'outputFolder')
    if ~exist(sPlotParams.outputFolder, 'dir')
        mkdir(sPlotParams.outputFolder)
    end
    simPrefix = strcat(sPlotParams.sDataset.actualDataDist, num2str(sPlotParams.sDataset.dim), ...
        'd', '_', sPlotParams.matrixForEigs);

    saveas(fig,strcat(sPlotParams.outputFolder, filesep, simPrefix, '_graphSigErr_', ...
        num2str(xInd(1)), '_to_', num2str(xInd(2))), 'epsc');
end
end