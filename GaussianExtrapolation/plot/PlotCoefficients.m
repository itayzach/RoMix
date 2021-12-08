function [] = PlotCoefficients(sSimParams, c, vLambda)

x0     = 10;
y0     = 250;
width  = 600;
height = 400;

cSquared = c.^2;

fig = figure;
tiledlayout('flow');
nexttile;
plot(cSquared, 'LineWidth', 2, 'DisplayName', '$c_m^2$');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14);
hold on;
plot(vLambda, 'LineWidth', 2, 'DisplayName', '$\lambda_m$');
plot(cSquared.^2 ./ vLambda, 'LineWidth', 2, 'DisplayName', '$c_m^2/\lambda_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'southwest')
% small_c_idx = find(abs(c_sq) < 1e-5);
% scatter(small_c_idx, c_sq(small_c_idx), 'filled');
set(gca,'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');

nexttile;
plot(c, 'LineWidth', 2, 'DisplayName', '$c_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'southwest')
set(gca,'FontSize', 14);

set(gcf,'Position', [x0 y0 width height])
if isfield(sSimParams, 'outputFolder')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_rkhs_norm_decay'], 'epsc');
end
end