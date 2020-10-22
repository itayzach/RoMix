function [] = PlotGraphToGraphTransform(G, G_tilde, f, f_tilde, G_title, G_tilde_title)

param.show_edges = false;
vPlotIndexes = (1:5)';(1:3:20)';
vNumAsText = strcat('\leftarrow ', {' '}, cellstr(num2str(vPlotIndexes)))';
v = G.coords;
G_dataDim = size(G.coords,2);
G_tilde_dataDim = size(G_tilde.coords,2);
v_tilde = G_tilde.coords;
figure;
subplot(1,2,1);
    if exist('f', 'var')
        gsp_plot_signal(G,f,param);
        title(['Graph-signal on $G$ ' newline G_title], 'Interpreter', 'latex', 'FontSize', 14);
    else
        gsp_plot_graph(G,param);
        title('Original graph', 'Interpreter', 'latex', 'FontSize', 14);
    end
    hold on;
    if G_dataDim == 3
        text(v(vPlotIndexes,1) + 0.05, v(vPlotIndexes,2)+0.01, v(vPlotIndexes,3), vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter3(v(vPlotIndexes,1), v(vPlotIndexes,2), v(vPlotIndexes,3), 'ko')
    else
        text(v(vPlotIndexes,1) + 0.05, v(vPlotIndexes,2)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v(vPlotIndexes,1), v(vPlotIndexes,2), 'ko')
    end
    view(0,90)

subplot(1,2,2);
    if exist('f_tilde', 'var')
        gsp_plot_signal(G_tilde,f_tilde,param);
        title(['Graph-signal on $\tilde{G}$' newline G_tilde_title], 'Interpreter', 'latex', 'FontSize', 14);
    else
        gsp_plot_graph(G_tilde,param);
        title('Transformed graph', 'Interpreter', 'latex', 'FontSize', 14);
    end
    hold on;
    if G_tilde_dataDim == 3
        text(v_tilde(vPlotIndexes,1) + 0.05, v_tilde(vPlotIndexes,2)+0.01, v_tilde(vPlotIndexes,3), vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter3(v_tilde(vPlotIndexes,1), v_tilde(vPlotIndexes,2), v_tilde(vPlotIndexes,3), 'ko')
    else
        text(v_tilde(vPlotIndexes,1) + 0.05, v_tilde(vPlotIndexes,2)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v_tilde(vPlotIndexes,1), v_tilde(vPlotIndexes,2), 'ko')
    end
set(gcf,'Position', [400 400 1000 400])  
end