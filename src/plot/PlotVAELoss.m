function PlotVAELoss(sPlotParams, history)
    sHistTmp = struct(history.history);
    fields = fieldnames(sHistTmp);
    nFeilds = numel(fields);
    for f = 1:nFeilds
        sHist.(fields{f}) = cellfun(@double,cell(sHistTmp.(fields{f})))';
    end
    sHist.epochs = cellfun(@double,cell(history.epoch))';

    windowStyle = get(0,'DefaultFigureWindowStyle');
    set(0,'DefaultFigureWindowStyle','normal')
    fig = figure('Name','Loss'); 
    hold on;
    plot(sHist.epochs, sHist.loss,'LineWidth', 2, 'DisplayName','Total');
    plot(sHist.epochs, sHist.kl_loss,'LineWidth', 2,'DisplayName','KL');
    plot(sHist.epochs, sHist.reconstruction_loss,'LineWidth', 2,'DisplayName','Reconstruction');

    legend('Location','northeast','Interpreter','latex','FontSize',14)
    xlabel('Epoch', 'Interpreter', 'latex', 'FontSize', 14)
    ylabel('Loss', 'Interpreter', 'latex', 'FontSize', 14)
    %% Size
    x0     = 10;
    y0     = 50;
    height = 350;
    width  = 600;
    set(gcf,'Position', [x0 y0 width height])
    %% Save
    figName = 'VAELoss';
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
    
    set(0,'DefaultFigureWindowStyle',windowStyle)
end