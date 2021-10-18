function PlotCoeffsMatrix(C, Cint)
figure('Name', 'C');
tiledlayout('flow')
nexttile;

imagesc(C);
colormap('jet');
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}_n^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
if exist('Cint', 'var')
    nexttile;
    imagesc(Cint);
    colormap('jet');
    colorbar();
    title('${\bf C^{{\bf int}}} = \tilde{{\bf \Phi}}_N^\dagger {\bf V^{{\bf int}}} $', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
end
end