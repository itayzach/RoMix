function [] = PlotPolyCdfDemonstration2(xMin, xMax, pCdf, invpCdf, muTilde, sigmaTilde)
nTestPoints = 50;
xSmallTest = linspace(xMin+0.5,xMax-0.5,nTestPoints)';

xTestTilde = T(pCdf, true, muTilde, sigmaTilde, xSmallTest);
xSmallTestEst = invT(invpCdf, muTilde, sigmaTilde, xTestTilde);

cmap = xSmallTest;
fig = figure('Name', 'Demonstrate T #2');
subplot(2,1,1)
    scatter(xSmallTest, zeros(1,nTestPoints), 100, cmap, 'o')
    hold on;
    scatter(xSmallTestEst, zeros(1,nTestPoints), 50, cmap, 'filled')
    colormap('jet')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    legend('$x_{{\bf test}}$', '$T^{-1}(T(x_{{\bf test}}))$','interpreter', 'latex', 'FontSize', 14);
    title('Original nodes','interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
subplot(2,1,2);
    scatter(xTestTilde, zeros(1,nTestPoints), 50, cmap, 'filled')
    colormap('jet')
    xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
    title('Transformed nodes','interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
% saveas(fig,strcat(outputFolder, filesep, 'fig4_x_to_xTilde'), figSaveType);
end