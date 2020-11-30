function [] = PlotGraphToGraphTransform(sSimParams, G1, G1_title, G2, G2_title, G3, G3_title, f1, f1_title, f2, f2_title, f3, f3_title, sampleInd)

nVerteciesText = 4;
param.show_edges = false;
if exist('sampleInd', 'var')
    vTextPlotIndexes = sampleInd(1:nVerteciesText);
    if ~sSimParams.b_plotSamplingPointsMarkers
        sampleInd = sampleInd(1:nVerteciesText);
    end
else
    rperm = randperm(size(G1.coords,1));
    vTextPlotIndexes = rperm(1:nVerteciesText)'; %(1:5)';(1:3:20)';
    sampleInd = vTextPlotIndexes;
end
vNumAsText = strcat('\leftarrow ', {' '}, cellstr(num2str(vTextPlotIndexes)))';

if exist('G3', 'var') && ~isempty(G3)
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
    if isempty(G)
        continue
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
        
        if sSimParams.b_plotSamplingPointsMarkers
            if G_dataDim == 3
                scatter3(v(sampleInd,1), v(sampleInd,2), v(sampleInd,3), 'ko')
                text(v(vTextPlotIndexes,1)*1.05, v(vTextPlotIndexes,2)*1.01, v(vTextPlotIndexes,3), vNumAsText,...
                    'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
            elseif G_dataDim == 2
                scatter(v(sampleInd,1), v(sampleInd,2), 50, 'ko')
                text(v(vTextPlotIndexes,1)*1.05, v(vTextPlotIndexes,2)*1.01, vNumAsText,...
                    'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
            elseif G_dataDim == 1
                shiftright = (max(v(:,1))-min(v(:,1)))/length(v);
                shiftup = (max(f)-min(f))/length(f);
                scatter(v(sampleInd), f(sampleInd), 'ko')
                text(v(vTextPlotIndexes)+2*shiftright, f(vTextPlotIndexes)+2*shiftup, vNumAsText,...
                    'FontWeight','bold','Color', 'black', 'BackgroundColor', 'white', 'Margin', 1,'EdgeColor','black')
            else
                error('invalid dim');
            end
        end
        view(0,90)
        N = size(v,1);
        plot_title = ['Graph-signal on ' G_title ' ($N$ = ' num2str(N) ')' newline f_title];
        if iGraph == nGraphs && length(f1) == length(f)
            err = norm(f1-f)/norm(f1);
            plot_title = strcat(plot_title, [ newline 'Error = ', num2str(err, '%d')]);
        end
        title(plot_title, 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
    else
        gsp_plot_graph(G,param);
        title(G_title, 'Interpreter', 'latex', 'FontSize', 14);
    end
end   
set(gcf,'Position', [200 400 nGraphs*500 400])

end