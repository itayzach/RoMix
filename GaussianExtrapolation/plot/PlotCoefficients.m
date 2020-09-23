function [] = PlotCoefficients(sSimParams, c, vLambda)

x0     = 10;
y0     = 250;
width  = 600;
height = 400;

fig = figure;
plot(vLambda, 'LineWidth', 2, 'DisplayName', '$\lambda_m$');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14);
hold on;
cSquared = c.^2;
plot(cSquared, 'LineWidth', 2, 'DisplayName', '$c_m^2$');
plot(cSquared.^2 ./ vLambda, 'LineWidth', 2, 'DisplayName', '$c_m^2/\lambda_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'southwest')
% small_c_idx = find(abs(c_sq) < 1e-5);
% scatter(small_c_idx, c_sq(small_c_idx), 'filled');
set(gca,'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');
set(gcf,'Position', [x0 y0 width height])
saveas(fig,[sSimParams.outputFolder filesep 'fig_rkhs_norm_decay'], 'epsc');
end