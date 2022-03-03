function PlotGraphMatrices(G, b_normlizedLaplacian)
figure; 
subplot(131)
    imagesc(G.W); 
    colorbar; 
    title('Kerel (weights) matrix', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(132); 
    imagesc(diag(G.d)); 
    colorbar; 
    title('Degree matrix', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);      
subplot(133); 
    imagesc(G.L); 
    colorbar; 
    if b_normlizedLaplacian
        title('Normalized Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
    else
        title('Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
    end
    set(gca,'FontSize', 14);    
set(gcf,'Position', [100 200 1800 400])
end

