function [] = PlotEigenDiffs(sSimParams, sDataset, nysRatio, vAnaVsNum, vNysVsNum)

%% eigen index
vM = 0:length(vAnaVsNum)-1;

%% Plot
fig = figure;
plot(vM, vAnaVsNum, 'LineWidth', 2, ...
    'DisplayName',  '$\sqrt{\frac{1}{T} \sum_{t=1}^T \big\| \phi_m^{Ana}(t) - \phi_m^{Num}(t) \big\|^2_2}$' );
hold on
plot(vM, vNysVsNum, 'LineWidth', 2, ...
    'DisplayName',  '$\sqrt{\frac{1}{T} \sum_{t=1}^T \big\| \phi_m^{Nys}(t) - \phi_m^{Num}(t) \big\|_2}$' );
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
title(strcat('$nNys = $',  num2str(nysRatio*sDataset.nTotal), '$\quad N = $', num2str(sDataset.nTotal), ...
    '$\quad \bigr(\frac{nNys}{N} = $', num2str(nysRatio, '%.2f'), '$\bigr)$'), ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
simPrefix = strcat(sDataset.actualDataDist, '_', num2str(nysRatio*100, '%d'), 'prec');
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_eigen_diffs_N_', num2str(sDataset.nTrain), '_nNys_', num2str(nysRatio*sDataset.nTrain)), 'epsc');


end