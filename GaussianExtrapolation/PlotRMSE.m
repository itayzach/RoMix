function [] = PlotRMSE(sSimParams, sDataset, nysRatio, vAnaVsNum, vNysVsNum)

%% eigen index
vM = 0:length(vAnaVsNum)-1;

%% Plot
fig = figure('Name', 'RMSE');
plot(vM, vAnaVsNum, 'LineWidth', 2, ...
    'DisplayName',  'RMSE$(\phi_m,  v_m)$' );
hold on
plot(vM, vNysVsNum, 'LineWidth', 2, ...
    'DisplayName',  'RMSE$(\hat{v}_m,  v_m)$' );
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
simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd', '_', num2str(nysRatio*100, '%d'), 'prec');
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_RMSE_N_', num2str(sDataset.nTrain), '_nNys_', num2str(nysRatio*sDataset.nTrain), '_', 'eigs_', num2str(length(vAnaVsNum))), 'epsc');


end