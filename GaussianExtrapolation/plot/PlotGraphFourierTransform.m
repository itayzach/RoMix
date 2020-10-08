function PlotGraphFourierTransform(sSimParams,G,f_hat)
param.show_edges = false;
figure; 
subplot(1,2,1)
    gsp_plot_signal_spectral(G,f_hat,param);
    xlabel('$\lambda_k$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\hat{f}$', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(1,2,2)
    stem(f_hat);
    xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$\hat{f}$', 'Interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
sgtitle('Graph signal spectrum', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf,'Position', [100 200 1200 400])
end

