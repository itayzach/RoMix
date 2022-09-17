function fig = PlotSpectrum(sPlotParams, vNysRatio, vLambda1, vLambda2, mLambda3, lambda1Str, lambda2Str, lambda3Str, figTitle)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

if ~exist('figTitle', 'var')
    figTitle = 'Eigenvalues';
end

fig = figure('Name', 'Spectrum');
M = length(vLambda1);
% subplot(2,1,1);
stem(0:M-1, vLambda1(1:M), 'filled', 'DisplayName', ['$' lambda1Str '$']);
if ~isempty(vLambda2)
    hold on;
    stem(0:M-1, vLambda2(1:M), 'DisplayName', ['$' lambda2Str '$']);
end
if ~isempty(vNysRatio)
    for r = 1:length(vNysRatio)
        nysRatio = vNysRatio(r);
        stem(0:M-1, mLambda3(r,:), '*', 'DisplayName', [['$' lambda3Str '$'] ' (' num2str(nysRatio*100, '%d') '\%)']);
    end
end
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0 M-1])
title(figTitle, 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

% subplot(2,1,2);
% stem(0:M-1, log10(abs(vLambda1(1:M))),'x', 'LineStyle','none', 'DisplayName', ['$\log|' lambda1Str '|$']);
% if ~isempty(vLambda2)
%     hold on;
%     stem(0:M-1, log10(abs(vLambda2(1:M))),'LineStyle','none', 'DisplayName', ['$\log|' lambda2Str '|$']);
% end
% if ~isempty(vNysRatio)
%     for r = 1:length(vNysRatio)
%         nysRatio = vNysRatio(r);
%         stem(0:M-1, mLambda3(r,:),'*', 'LineStyle','none', 'DisplayName', [['$\log(' lambda3Str ')$'] ' (' num2str(nysRatio*100, '%d') '\%)']);
%     end
% end
set(gca, 'YScale', 'log')
% legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
% xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% xlim([0 M + 8])
% title([figTitle, ' (log-scale)'], 'Interpreter', 'latex', 'FontSize', 14);
% set(gca,'FontSize', 14);

set(gcf,'Position',[-1000 100 600 200])
set(0,'DefaultFigureWindowStyle',windowStyle)
%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'spectrum';
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end