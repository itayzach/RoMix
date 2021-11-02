function PlotGraphSignalAnalysis(sig, sigRecPhi, sigRecV, sigRef, sigInt, sigNys, sigHatPhi, sigHatV, sigRefHatPhi)

n = length(sig);
N = length(sigRef);

figure('Name', 'Graph signals analysis');
tiledlayout('flow')

nexttile;
stem(sigHatPhi, 'filled', 'DisplayName', '${{\bf \Phi}}_n^\dagger s$');
hold on
stem(sigRefHatPhi, 'filled', 'DisplayName', '${{\bf \Phi}}_N^\dagger s^{{\bf ref}}$');
%     stem(C*sigHatV, 'filled', 'DisplayName', '${\bf C} {\bf V}^T s$');
%     stem(sigHatPhi2, 'filled', 'DisplayName', '$?$');
title('$\Phi$ coeffs', 'interpreter', 'latex', 'FontSize', 14)
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
set(gca,'FontSize', 14);

nexttile;
stem(sigHatV, 'filled');
title('$V$ coeffs', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(1:n, abs(sig-sigRecPhi), 'x', 'DisplayName', '$|s-s_{\Phi}^{{\bf rec}}|$');
hold on;
plot(1:n, abs(sig-sigRecV), '+', 'DisplayName', '$|s-s_{V}^{{\bf rec}}|$');
plot(1:n, mean(abs(sig-sigRecPhi))*ones(1,n), 'LineWidth', 3, 'DisplayName', '${{\bf avg}}({|s-s_{\Phi}^{{\bf rec}}|})$');
plot(1:n, mean(abs(sig-sigRecV))*ones(1,n), 'LineWidth', 3, 'DisplayName', '${{\bf avg}}({|s-s_{V}^{{\bf rec}}|})$');
xlim([1,n]);
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title(['Projection error', newline, '($n$ given nodes)'], 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(n+1:N, abs(sigRef(n+1:end)-sigInt(n+1:end)), 'x', 'DisplayName', '$|s^{{\bf ref}}-s^{{\bf int}}|$');
hold on;
plot(n+1:N, abs(sigRef(n+1:end)-sigNys(n+1:end)), '+', 'DisplayName', '$|s^{{\bf ref}}-s_{{\bf nys}}|$');
plot(n+1:N, mean(abs(sigRef(n+1:end)-sigInt(n+1:end)))*ones(1,N-n), 'LineWidth', 3, 'DisplayName', '${{\bf avg}}({|s^{{\bf ref}}-s^{{\bf int}}|})$');
plot(n+1:N, mean(abs(sigRef(n+1:end)-sigNys(n+1:end)))*ones(1,N-n), 'LineWidth', 3, 'DisplayName', '${{\bf avg}}({|s^{{\bf ref}}-s_{{\bf nys}}|})$');
if n+1 < N
    xlim([n+1,N]);
end
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title(['Interpolation error', newline, '($N-n$ nodes)'], 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(1:n, sig, 'DisplayName', '$s$');
hold on;
plot(1:n, sigRecPhi, 'DisplayName', '$s_{\Phi}^{{\bf rec}}$');
%     plot(1:n, sigRecV,  'DisplayName', '$s_{V}^{{\bf rec}}$');
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title('graph signals on $n$ nodes', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(1:N, sigRef, 'DisplayName', '$s^{{\bf ref}}$');
hold on;
plot(1:N, sigInt, 'DisplayName', '$s^{{\bf int}}$');
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title('graph signals on $N$ nodes', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

nexttile;
plot(1:N, sigRef, 'DisplayName', '$s^{{\bf ref}}$');
hold on;
plot(1:N, sigNys, 'DisplayName', '$s^{{\bf nys}}$');
legend('interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',2)
title('graph signals on $N$ nodes', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

x0     = 10;
y0     = 50;
height = 800;
width  = 1500;
set(gcf,'Position', [x0 y0 width height])
end