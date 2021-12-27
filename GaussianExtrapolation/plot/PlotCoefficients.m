function [] = PlotCoefficients(sSimParams, c, vLambda, cRep)

x0     = 10;
y0     = 250;
width  = 600;
height = 400;

fig = figure;
tiledlayout('flow');
nexttile;
plot(abs(c), 'LineWidth', 2, 'DisplayName', '$|c_m|$');
plot(c.^2, 'LineWidth', 2, 'DisplayName', '$c_m^2$');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14);
hold on;
plot(vLambda, 'LineWidth', 2, 'DisplayName', '$\lambda_m$');
% plot(c.^2 ./ vLambda, 'LineWidth', 2, 'DisplayName', '$c_m^2/\lambda_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthOutside', 'NumColumns', 3)
% small_c_idx = find(abs(c_sq) < 1e-5);
% scatter(small_c_idx, c_sq(small_c_idx), 'filled');
set(gca,'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');

nexttile;
plot(c, 'LineWidth', 2, 'DisplayName', '$c_m$');
if exist('cRep', 'var') && ~isempty(cRep)
    hold on;
    plot(cRep, 'LineWidth', 2, 'DisplayName', '$\hat{c}_m$');
end
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthOutside', 'NumColumns', 2)
set(gca,'FontSize', 14);

if exist('cRep', 'var') && ~isempty(cRep)
    nexttile;
    plot(abs(c-cRep), 'LineWidth', 2, 'DisplayName', '$|c_m-\hat{c}_m|$');
    legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'SouthOutside', 'NumColumns', 2)
    set(gca,'FontSize', 14);
end

set(gcf,'Position', [x0 y0 width height])
if isfield(sSimParams, 'outputFolder')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_rkhs_norm_decay'], 'epsc');
end
end