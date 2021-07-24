function [] = PlotRMSE(sPlotParams, vRmseInt, vRmseNys)

M = length(vRmseInt);

figure('Name', 'RMSE');
plot((0:M-1)', vRmseInt.', '-o', 'LineWidth', 2, 'DisplayName', ...
    'RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$'), hold on;
plot((0:M-1)', vRmseNys.', '-x', 'LineWidth', 2, 'DisplayName', ...
    'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$')
rmseylim = max([vRmseInt vRmseNys]);
ylim([0 rmseylim]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);

%% Save
if isfield(sPlotParams, 'outputFolder')
    if ~exist(sPlotParams.outputFolder, 'dir')
        mkdir(sPlotParams.outputFolder)
    end
    saveas(fig,strcat(sPlotParams.outputFolder, filesep, 'RMSE_', 'eigs_', num2str(length(vRmseInt))), 'epsc');
end


end