function PlotVerticesTransformation(sSimParams,dataDim,v,v_tilde,graphName,verticesTransformation)
if dataDim == 1
    histPlotRows = 1;
else
    histPlotRows = 3;
end

fig = figure; 
subplot(histPlotRows,2,1)
histogram(v(:,1),100);
set(gca,'FontSize', 14);
title('$v_1$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
subplot(histPlotRows,2,2)
histfit(v_tilde(:,1),100);
set(gca,'FontSize', 14);
title('$\tilde{v}_1$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
if dataDim == 2 || dataDim == 3
    subplot(histPlotRows,2,3)
    histogram(v(:,2),100);
    title('$v_2$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca,'FontSize', 14);
    subplot(histPlotRows,2,4)
    histfit(v_tilde(:,2),100);
    title('$\tilde{v}_2$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca,'FontSize', 14);
    
    if dataDim == 2
        subplot(histPlotRows,2,5)
        hist3(v,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
        colormap(gca, 'hot')
        colorbar()
        view(2)
        xlim([ min(v(:,1)) max(v(:,1))])
        ylim([ min(v(:,2)) max(v(:,2))])
        title('$v$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);

        subplot(histPlotRows,2,6)
        hist3(v_tilde,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
        colormap(gca, 'hot')
        colorbar()
        view(2)
        xlim([ min(v_tilde(:,1)) max(v_tilde(:,1))])
        ylim([ min(v_tilde(:,2)) max(v_tilde(:,2))])
        title('$\tilde{v}$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
    else
        subplot(histPlotRows,2,5)
        histogram(v(:,3),100);
        title('$v_3$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
        subplot(histPlotRows,2,6)
        histfit(v_tilde(:,3),100);
        title('$\tilde{v}_3$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
    end    
end
set(gcf,'Position', [100 200 1000 histPlotRows*250])
if isfield(sSimParams, 'outputFolder')
    saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_hists_transformed'), 'epsc');
end
end

