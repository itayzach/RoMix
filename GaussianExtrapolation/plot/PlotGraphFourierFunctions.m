function PlotGraphFourierFunctions(G, nEigs)
param.show_edges = false;
figure;
for i=1:min(10,nEigs)
    subplot(2,5,i)
    gsp_plot_signal(G,G.U(:,i),param); 
    view(0,90)
    title(['u_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(G.e(i), '%.4f')]);
end
sgtitle('Laplacian eigenvectors');
set(gcf,'Position', [200 200 1600 600])
end
