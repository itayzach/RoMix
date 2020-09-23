function fig = PlotSpectrum(sSimParams, sDataset, vNysRatio, vLambdaAnalytic, vLambdaNumeric, mLambdaNystrom)

fig = figure('Name', 'Spectrum');
M = length(vLambdaAnalytic);
subplot(2,1,1);
stem(0:M-1, vLambdaAnalytic, 'x', 'DisplayName', '$\lambda^{\phi}_m$');
hold on;
stem(0:M-1, vLambdaNumeric, 'DisplayName', '$\lambda^{v}_m$');
if ~isempty(vNysRatio)
    for r = 1:length(vNysRatio)
        nysRatio = vNysRatio(r);
        stem(0:M-1, mLambdaNystrom(r,:), '*', 'DisplayName', ['$\lambda^{\hat{v}}_m$' ' (' num2str(nysRatio*100, '%d') '\%)']);
    end
end
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0 M + 8])
title('Eigenvalues', 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

subplot(2,1,2);
stem(0:M-1, log(vLambdaAnalytic),'x', 'LineStyle','none', 'DisplayName', '$\log(\lambda^{\phi}_m)$');
hold on;
stem(0:M-1, log(vLambdaNumeric),'LineStyle','none', 'DisplayName', '$\log(\lambda^{v}_m)$');
if ~isempty(vNysRatio)
    for r = 1:length(vNysRatio)
        nysRatio = vNysRatio(r);
        stem(0:M-1, log(mLambdaNystrom(r,:)),'*', 'LineStyle','none', 'DisplayName', ['$\log(\lambda^{\hat{v}}_m)$' ' (' num2str(nysRatio*100, '%d') '\%)']);
    end
end

legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0 M + 8])
title('Eigenvalues log-scale', ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
set(gcf,'Position',[100 100 600 500])

%% Save
try
    if ~exist(sSimParams.outputFolder, 'dir')
        mkdir(sSimParams.outputFolder)
    end

    simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd');
    saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_spectrum'), 'epsc');
catch
end
end