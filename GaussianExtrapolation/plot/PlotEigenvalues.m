function fig = PlotEigenvalues(sSimParams, actualDataDist, lambda1_title, vLambda1, lambda2_title, vLambda2)

fig = figure('Name', 'Spectrum');
M1 = length(vLambda1);
M2 = length(vLambda2);
subplot(2,1,1);
plot(0:M1-1, vLambda1(1:M1), 'x', 'DisplayName', lambda1_title);
if exist('vLambda2','var') && ~isempty(vLambda2)
    hold on;
    plot(0:M2-1, vLambda2(1:M2), 'o', 'DisplayName', lambda2_title);
end

xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
% ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0 max(M1,M2) + 1])
title('Eigenvalues', 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

subplot(2,1,2);
semilogy(0:M1-1, vLambda1(1:M1),'x', 'LineStyle','none', 'DisplayName', lambda1_title);
if exist('vLambda2','var') && ~isempty(vLambda2)
    hold on;
    semilogy(0:M2-1, vLambda2(1:M2),'o','LineStyle','none', 'DisplayName', lambda2_title);
end

legend('Interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
xlim([0 max(M1,M2) + 1])
title('Eigenvalues (log-scale)', 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
set(gcf,'Position',[600 100 600 800])

%% Save
if isfield(sSimParams, 'outputFolder')
    if ~exist(sSimParams.outputFolder, 'dir')
        mkdir(sSimParams.outputFolder)
    end
    saveas(fig,strcat(sSimParams.outputFolder, filesep, actualDataDist, '_spectrum'), 'epsc');
end
end