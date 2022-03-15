function PlotGraphSignalAnalysis(sig, sigRecPhi, sigRecV, sigRef, sigInt, sigNys, sigHatPhi, sigHatV, sigRefHatPhi)

n = length(sig);
N = length(sigRef);

figure('Name', 'Graph signals analysis');
tiledlayout('flow')

nexttile;
stem(sigHatPhi, 'filled', 'DisplayName', '${{\bf \Phi}}_n^\dagger s$');
hold on
stem(sigRefHatPhi, 'filled', 'DisplayName', '${{\bf \Phi}}_N^\dagger s^{{\bf gt}}$');
title('$\Phi$ coeffs', 'interpreter', 'latex', 'FontSize', 14)
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
set(gca,'FontSize', 14);

nexttile;
stem(sigHatV, 'filled');
title('$V$ coeffs', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);


nexttile;
plot(1:n, sig, 'DisplayName', '$s$');
hold on;
plot(1:n, sigRecPhi, 'o', 'DisplayName', '$s_{\Phi}^{{\bf rec}}$');
plot(1:n, sigRecV,  'o', 'DisplayName', '$s_{V}^{{\bf rec}}$');
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title('graph signals on $n$ nodes', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(1:N, sigRef, 'DisplayName', '$\tilde{s}$');
hold on;
plot(1:N, sigInt, 'o', 'DisplayName', '$\tilde{s}^{{\bf RoMix}}$');
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title('graph signals on $N$ nodes', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(1:N, sigRef, 'DisplayName', '$\tilde{s}$');
hold on;
plot(1:N, sigNys, 'o', 'DisplayName', '$\tilde{s}^{{\bf nys}}$');
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title('graph signals on $N$ nodes', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

x0     = 10;
y0     = 50;
height = 800;
width  = 1500;
set(gcf,'Position', [x0 y0 width height])
end