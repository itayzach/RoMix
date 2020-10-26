function [] = PlotGraphToGraphTransformAndBack(G1, G2, G3, G1_title, G2_title, G3_title, f1, f2, f3, f1_title, f2_title, f3_title)

param.show_edges = false;
vPlotIndexes = (1:5)';(1:3:20)';
vNumAsText = strcat('\leftarrow ', {' '}, cellstr(num2str(vPlotIndexes)))';
v = G1.coords;
G_dataDim = size(G1.coords,2);
G_tilde_dataDim = size(G2.coords,2);
v_tilde = G2.coords;
figure;
subplot(1,3,1);
    if exist('f1', 'var')
        if size(v,2) > 1
            gsp_plot_signal(G1,f1,param);
        else
            plot(G1.coords, f1, '.');
        end
        title(['Graph-signal on ' G1_title newline f1_title], 'Interpreter', 'latex', 'FontSize', 14);
    else
        gsp_plot_graph(G1,param);
        title('Original graph', 'Interpreter', 'latex', 'FontSize', 14);
    end
    hold on;
    if G_dataDim == 3
        text(v(vPlotIndexes,1) + 0.05, v(vPlotIndexes,2)+0.01, v(vPlotIndexes,3), vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter3(v(vPlotIndexes,1), v(vPlotIndexes,2), v(vPlotIndexes,3), 'ko')
    elseif G_dataDim == 2
        text(v(vPlotIndexes,1) + 0.05, v(vPlotIndexes,2)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v(vPlotIndexes,1), v(vPlotIndexes,2), 'ko')
    elseif G_dataDim == 1
        text(v(vPlotIndexes) + 0.05, f1(vPlotIndexes)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v(vPlotIndexes), f1(vPlotIndexes), 'ko')
    else
        error('invalid dim');
    end
    view(0,90)

subplot(1,3,2);
    if exist('f2', 'var')
        if size(v,2) > 1
            gsp_plot_signal(G2,f2,param);
        else
            plot(G2.coords, f2, '.');
        end
        title(['Graph-signal on ' G2_title  newline f2_title], 'Interpreter', 'latex', 'FontSize', 14);
    else
        gsp_plot_graph(G2,param);
        title('Transformed graph', 'Interpreter', 'latex', 'FontSize', 14);
    end
    hold on;
    if G_tilde_dataDim == 3
        text(v_tilde(vPlotIndexes,1) + 0.05, v_tilde(vPlotIndexes,2)+0.01, v_tilde(vPlotIndexes,3), vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter3(v_tilde(vPlotIndexes,1), v_tilde(vPlotIndexes,2), v_tilde(vPlotIndexes,3), 'ko')
    elseif G_dataDim == 2
        text(v_tilde(vPlotIndexes,1) + 0.05, v_tilde(vPlotIndexes,2)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v_tilde(vPlotIndexes,1), v_tilde(vPlotIndexes,2), 'ko')
    elseif G_dataDim == 1
        text(v_tilde(vPlotIndexes,1) + 0.05, f2(vPlotIndexes)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v_tilde(vPlotIndexes), f2(vPlotIndexes), 'ko')
    else
        error('invalid dim');
    end
subplot(1,3,3);
    if exist('f3', 'var')
        if size(v,2) > 1
            gsp_plot_signal(G3,f3,param);
        else
            plot(G3.coords, f3, '.');
        end
        err = norm(f1-f3)/norm(f1);
        title(['Graph-signal on ' G3_title  newline f3_title newline 'error = ' num2str(err)], 'Interpreter', 'latex', 'FontSize', 14);
    else
        gsp_plot_graph(G3,param);
        title('Transformed graph', 'Interpreter', 'latex', 'FontSize', 14);
    end
    hold on;
    if G_dataDim == 3
        text(v(vPlotIndexes,1) + 0.05, v(vPlotIndexes,2)+0.01, v(vPlotIndexes,3), vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter3(v(vPlotIndexes,1), v(vPlotIndexes,2), v(vPlotIndexes,3), 'ko')
    elseif G_dataDim == 2
        text(v(vPlotIndexes,1) + 0.05, v(vPlotIndexes,2)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v(vPlotIndexes,1), v(vPlotIndexes,2), 'ko')
    elseif G_dataDim == 1
        text(v(vPlotIndexes) + 0.05, f3(vPlotIndexes)+0.01, vNumAsText,...
            'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
        scatter(v(vPlotIndexes), f3(vPlotIndexes), 'ko')
    else
        error('invalid dim');
    end
    view(0,90)
    
set(gcf,'Position', [200 400 1500 400])  
end