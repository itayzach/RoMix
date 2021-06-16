function [] = PlotPolyCdfDemonstration2(xMin, xMax, pCdf, invpCdf, muTilde, sigmaTilde, xTrain, b_kde, kde_bw)
if ~exist('b_kde', 'var')
    b_kde = false;
end
if ~exist('kde_bw', 'var')
    kde_bw = [];
end

nTestPoints = 50;
xSmallTest = linspace(xMin+0.5,xMax-0.5,nTestPoints)';
xTestTilde = T(pCdf, true, muTilde, sigmaTilde, xSmallTest, xTrain, b_kde, kde_bw);
if b_kde
else
    xSmallTestEst = Tinv(invpCdf, muTilde, sigmaTilde, xTestTilde);
end


cmap = xSmallTest;
fig = figure('Name', 'Demonstrate T #2');
subplot(2,1,1)
    if b_kde 
        scatter(xSmallTest, zeros(1,nTestPoints), 50, cmap, 'filled', 'DisplayName', '$x_{{\bf test}}$')
    else
        scatter(xSmallTest, zeros(1,nTestPoints), 100, cmap, 'o', 'DisplayName', '$x_{{\bf test}}$')
        hold on;
        scatter(xSmallTestEst, zeros(1,nTestPoints), 50, cmap, 'filled', 'DisplayName', '$T^{-1}(T(x_{{\bf test}}))$')
    end
    colormap('jet')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    legend('interpreter', 'latex', 'FontSize', 14);
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