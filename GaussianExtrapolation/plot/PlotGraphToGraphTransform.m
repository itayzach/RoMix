function [] = PlotGraphToGraphTransform(G1, G1_title, G2, G2_title, G3, G3_title, f1, f1_title, f2, f2_title, f3, f3_title)

param.show_edges = false;

rperm = randperm(size(G1.coords,1));
vPlotIndexes = rperm(1:5)'; %(1:5)';(1:3:20)';
vNumAsText = strcat('\leftarrow ', {' '}, cellstr(num2str(vPlotIndexes)))';

if exist('G3', 'var')
    nGraphs = 3;
else
    nGraphs = 2;
end

figure;
for iGraph = 1:nGraphs
    if iGraph == 1
        G = G1;
        G_title = G1_title;
        if exist('f1', 'var')
            f = f1;
            f_title = f1_title;
        end
        
    elseif iGraph == 2
        G = G2;
        G_title = G2_title;
        if exist('f2', 'var')
            f = f2;
            f_title = f2_title;
        end
    elseif iGraph == 3
        G = G3;
        G_title = G3_title;
        if exist('f3', 'var')
            f = f3;
            f_title = f3_title;
        end
    end
    v = G.coords;
    G_dataDim = size(G.coords,2);
    
    subplot(1,nGraphs,iGraph);
    if exist(['f' num2str(iGraph)], 'var')
        if size(v,2) > 1
            gsp_plot_signal(G,f,param);
        else
            plot(G.coords, f, '.');
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
            text(v(vPlotIndexes) + 0.05, f(vPlotIndexes)+0.01, vNumAsText,...
                'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
            scatter(v(vPlotIndexes), f(vPlotIndexes), 'ko')
        else
            error('invalid dim');
        end
        view(0,90)
        plot_title = ['Graph-signal on ' G_title newline f_title];
        if iGraph == 3
            err = norm(f1-f3)/norm(f1);
            plot_title = strcat(plot_title, [ newline 'Error = ', num2str(err, '%d')]);
        end
        title(plot_title, 'Interpreter', 'latex', 'FontSize', 14);
    else
        gsp_plot_graph(G,param);
        title(G_title, 'Interpreter', 'latex', 'FontSize', 14);
    end
end   
set(gcf,'Position', [200 400 nGraphs*500 400])

end