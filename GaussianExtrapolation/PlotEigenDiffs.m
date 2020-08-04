function [] = PlotEigenDiffs(sParams, vAnaVsNum, vNysVsNum)

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
title(strcat('$nNys = $',  num2str(sParams.nPointsNystrom), '$\quad N = $', num2str(sParams.nRandPoints), ...
    '$\quad \bigr(\frac{nNys}{N} = $', num2str(sParams.nysRatio, '%.2f'), '$\bigr)$'), ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

%% Save
if sParams.sSim.twomoons_dataset
    simPrefix = strcat('Two_moons_', num2str(sParams.nysRatio*100, '%d'), 'prec');
else
    simPrefix = strcat('Gaussian_', num2str(sParams.nysRatio*100, '%d'), 'prec');
end
if ~exist(sParams.sSim.outputFolder, 'dir')
    mkdir(sParams.sSim.outputFolder)
end

saveas(fig,strcat(sParams.sSim.outputFolder, filesep, simPrefix, '_eigen_diffs_N_', num2str(sParams.nRandPoints), '_nNys_', num2str(sParams.nPointsNystrom)), 'epsc');


end