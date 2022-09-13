function fig = PlotDataset(sPlotParams, sPreset, x, y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, b_transform, cXAxisLabels)
prevWindowStyle = get(0,'DefaultFigureWindowStyle');
if ~exist('windowStyle', 'var')
    windowStyle = prevWindowStyle;
end
set(0,'DefaultFigureWindowStyle',windowStyle)

if ~exist('b_transform', 'var')
    b_transform = false;
end

[n, dim] = size(x);
if exist('sPlotParams', 'var') && ~isempty(sPlotParams)
    actualDataDist = sPlotParams.actualDataDist;
else
    actualDataDist = '';
end

if ~exist('cXAxisLabels', 'var') || isempty(cXAxisLabels)
    cXAxisLabels{1} = '$z_1$';
    cXAxisLabels{2} = '$z_2$';
end

if ismember(actualDataDist, {'USPS', 'MNIST'})
    %if sDistParams.dim < 28*28 && strcmp(actualDataDist,'MNIST')
    if isfield(sDistParams, 'vae')
        % Transform from z to x
        x_train_rec = pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=x,vae=sDistParams.vae);
        xTrainRec = double(x_train_rec);
        plotInd = randperm(n,nGmmPoints);
        xDigits = reshape(xTrainRec(plotInd,:,:),[],28*28);
    else
        xDigits = x(randperm(n,nGmmPoints),:);
    end
    b_transpose = strcmp(sPlotParams.actualDataDist, 'MNIST') && sDistParams.dim == 28*28;
    PlotDigits(sPlotParams, xDigits, [], b_transpose, strcat(pltTitle, " (", actualDataDist, ")"), 'data');
    
    if exist('sDistParams', 'var') 
        [xGmm,~] = random(sDistParams.GMModel, nGmmPoints);
        if sDistParams.dim < 28*28 && strcmp(actualDataDist,'MNIST')
            % Transform from z to x
            x_gmm = pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=xGmm,vae=sDistParams.vae);
            xGmm = double(x_gmm);
            xGmm = reshape(xGmm,[],28*28);
        end
        PlotDigits(sPlotParams, xGmm, [], b_transpose, plt2Title, ['gmm_' num2str(sDistParams.estNumComponents) 'comp']);
    end
    
end

%% GMM
fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
compIdx(:,1) = ones(nGmmPoints,1);
if exist('sDistParams', 'var')
    b_plotDistModel = true;
    b_spectclust = isfield(sDistParams,'SCcompIdx');
    tiledlayout(1,2 + b_spectclust);
    vAx(2) = nexttile(2);
    [xGmm,compIdx(:,2)] = random(sDistParams.GMModel, nGmmPoints);
    if dim > 3 % Apply PCA
        pcaDim = 2;
        dim = pcaDim;
        [~,zGmm] = pca(xGmm);
        xGmm = zGmm(:,1:pcaDim);
        [~,z] = pca(x);
        x = z(:,1:pcaDim);
    end
    if b_transform && ismember(actualDataDist, {'SwissRoll'})
        sN = xGmm(:,1);
        a = sPreset.sDatasetParams.a;
        theta = linspace(0, sPreset.sDatasetParams.maxTheta, 100);
        s = SwissRollArclength(theta, a);
        thetaN = interp1(s, theta, sN);

        x1 = a * cos(thetaN) .* thetaN;
        x2 = a * sin(thetaN) .* thetaN;
        x3 = xGmm(:,2);
        xGmm = [x1, x2, x3];
    end
    xMax = max(x);
    xMin = min(x);
    for i=1:1+b_spectclust
        if i == 2
            xGmm = x;
            compIdx(:,3) = sDistParams.SCcompIdx;
            vAx(3) = nexttile(3);
        end
        if dim == 1
            %scatter(xGmm, zeros(1,nGmmPoints), 50, compIdx(:,i+1), 'filled')
            scatter(xGmm, zeros(1,nGmmPoints), 50, 'filled')
            xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
            set(gca,'YTick',[],'FontSize', 14);
            xlim([xMin(1), xMax(1)])
        elseif dim == 2
            if ismember(actualDataDist, {'SwissRoll'})
                %gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(sDistParams.GMModel,[x0 y0]),x,y);
                %xMax = max(cell2mat(sDistParams.mu')) + 3*max(cell2mat(sDistParams.sigma'));
                %xMin = min(cell2mat(sDistParams.mu')) - 3*min(cell2mat(sDistParams.sigma'));
                %fcontour(gmPDF, [xMin(1), xMax(1), xMin(2), xMax(2)]);
                %hold on;
                scatter3(xGmm(:,1),xGmm(:,2),compIdx(:,i+1),[],compIdx(:,i+1),'filled');
                colormap(lines(sDistParams.GMModel.NumComponents));
                h = colorbar('TickLabelInterpreter', 'latex');
                h.Limits = [0 sDistParams.GMModel.NumComponents];
                h.Ticks = (0.5:sDistParams.GMModel.NumComponents);
                h.TickLabels = 1:sDistParams.GMModel.NumComponents;
            else
                scatter(xGmm(:,1),xGmm(:,2),[],'filled');
            end
            grid on;
            if ismember(actualDataDist, {'SwissRoll'})
                xlabel(cXAxisLabels{1}, 'interpreter', 'latex', 'FontSize', 16);
                ylabel(cXAxisLabels{2}, 'interpreter', 'latex', 'FontSize', 16);
            else
                xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
                ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
            end
            view(2);
            set(gca,'FontSize', 14);
            xlim([xMin(1), xMax(1)])
            ylim([xMin(2), xMax(2)])
        elseif dim == 3
            if ismember(actualDataDist, {'SwissRoll'}) && sPreset.sDatasetParams.b_randn
                scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[],compIdx(:,i+1), 'filled');
                colormap(lines(sDistParams.GMModel.NumComponents));
                h = colorbar('TickLabelInterpreter', 'latex');
                h.Ticks = (1:sDistParams.GMModel.NumComponents);
                h.TickLabels = 1:sDistParams.GMModel.NumComponents;
            else
                scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[], 'filled');
            end
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
            zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
            view(30,75);
            set(gca,'FontSize', 14);
            xlim([xMin(1), xMax(1)])
            ylim([xMin(2), xMax(2)])
            zlim([xMin(3), xMax(3)])
        end
        if ~isempty(plt2Title)
            title(plt2Title, 'Interpreter', 'latex', 'FontSize', 14)
        end
    end
    vAx(1) = nexttile(1);
    UpdateCursorDataTip(fig, vAx, compIdx);
    
else
    xMax = max(x);
    xMin = min(x);
end

%% Dataset
if ismember(actualDataDist, {'Uniform', 'SwissRoll'}) && dim == 2
    hRatio = 2*(xMax(2)-xMin(2))/(xMax(1)-xMin(1));
else
    hRatio = 1;
end
if dim == 1
    scatter(x, zeros(1,n), 50, ones(1,n), 'filled')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
    xlim([xMin(1), xMax(1)])
elseif dim == 2
    if ~isempty(y) && isequal(y, floor(y)) && sum(y(:) > 0)
        y = ConvertSignalByDataset(actualDataDist, y)-1;
        scatter3(x(:,1), x(:,2), y, [], y, 'filled');
        if exist('ax', 'var')
            colormap(ax(2), jet(length(unique(y)))); 
        else
            colormap(jet(length(unique(y))));
        end
        colorbar;
    else
        scatter(x(:,1), x(:,2), 'filled');
        grid on;
    end
    if ismember(actualDataDist, {'SwissRoll'})
        xlabel(cXAxisLabels{1}, 'interpreter', 'latex', 'FontSize', 16);
        ylabel(cXAxisLabels{2}, 'interpreter', 'latex', 'FontSize', 16);
    else
        xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    end
    view(2);
    set(gca,'FontSize', 14);
    if exist('xMin', 'var') && exist('xMax', 'var')
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
    end
elseif dim == 3
    scatter3(x(:,1), x(:,2), x(:,3), 'filled');
    xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
%     view(10,5);
    view(30,75);
    set(gca,'FontSize', 14);
    if exist('xMin', 'var') && exist('xMax', 'var')
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
        zlim([xMin(3), xMax(3)])
    end
    hold on;
%     scatter3(sDistParams.GMModel.mu(:,1), sDistParams.GMModel.mu(:,2), sDistParams.GMModel.mu(:,3), 100, 'r', 'filled');
%     for c = 1:sDistParams.GMModel.NumComponents 
%         quiver3(sDistParams.u{c}(:,1), sDistParams.u{c}(:,2), sDistParams.u{c}(:,3),zeros(3,1),zeros(3,1),zeros(3,1));
%         hold on;
%     end
end
if ~isempty(pltTitle)
    title(strcat(pltTitle, " (", actualDataDist, ")"), 'Interpreter', 'latex', 'FontSize', 14)
end

%% Size
if strcmp(windowStyle, 'normal')
    x0     = 400;
    y0     = 400;
    height = hRatio*400;
    width  = 600*(1 + (exist('b_plotDistModel', 'var') && b_plotDistModel));
    set(gcf,'Position', [x0 y0 width height])
end
set(0,'DefaultFigureWindowStyle',prevWindowStyle)
%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder') && exist('fig','var')
    if exist('b_plotDistModel', 'var') && b_plotDistModel
        figName = 'dataset_vs_dist';
    else
        figName = 'dataset';
    end
    if ismember(actualDataDist, {'SwissRoll'})
        if dim == 3
            figName = [figName, '_x'];
        else
            figName = [figName, '_z'];
        end
    end
    set(fig,'renderer','Painters')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end