function [] = PlotRMSE(sSimParams, sDataset, vNysRatio, vAnaVsNum, mNysVsNum)

%% eigen index
vM = 0:length(vAnaVsNum)-1;

%% Plot
fig = figure('Name', 'RMSE');
plot(vM, vAnaVsNum, 'LineWidth', 2, 'DisplayName',  'RMSE$(\phi_m, v_m)$' );
hold on
for r = 1:length(vNysRatio)
    nysRatio = vNysRatio(r);
    plot(vM, mNysVsNum(r,:), 'LineWidth', 2, 'DisplayName',  ['RMSE$(\hat{v}_m, v_m)$' ' (' num2str(nysRatio*100, '%d') '\%)']);
end
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
% title(strcat('$nNys = $',  num2str(nysRatio*sDataset.nTotal), '$\quad N = $', num2str(sDataset.nTotal), ...
%     '$\quad \bigr(\frac{nNys}{N} = $', num2str(nysRatio, '%.2f'), '$\bigr)$'), ...
%     'Interpreter', 'latex', 'FontSize', 14);
ylim([0 1.5*max(max(mNysVsNum))])
set(gca,'FontSize', 14);

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end

simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd');
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_RMSE_', 'eigs_', num2str(length(vAnaVsNum))), 'epsc');


end