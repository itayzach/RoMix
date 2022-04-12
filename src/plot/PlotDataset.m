function fig = PlotDataset(sPlotParams, x, y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle)
prevWindowStyle = get(0,'DefaultFigureWindowStyle');
if ~exist('windowStyle', 'var')
    windowStyle = prevWindowStyle;
end
set(0,'DefaultFigureWindowStyle',windowStyle)

[n, dim] = size(x);
if exist('sPlotParams', 'var') && ~isempty(sPlotParams)
    actualDataDist = sPlotParams.actualDataDist;
else
    actualDataDist = '';
end

if ismember(actualDataDist, {'USPS', 'MNIST'})
    if sDistParams.dim < 28*28 && strcmp(actualDataDist,'MNIST')
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
        PlotDigits(sPlotParams, xGmm, [], b_transpose, plt2Title, 'gmm');
    end
    
end
if dim <= 3
    fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
    compIdx(:,1) = ones(n,1);
    %% GMM
    if exist('sDistParams', 'var')
        nGmmPoints = n;
        b_plotDistModel = true;
        b_spectclust = isfield(sDistParams,'SCcompIdx');
        tiledlayout(1,2 + b_spectclust);
        vAx(2) = nexttile(2);
        [xGmm,compIdx(:,2)] = random(sDistParams.GMModel, nGmmPoints);
        xMax = max([max(xGmm); max(x)]);
        xMin = min([min(xGmm); min(x)]);
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
                scatter(xGmm(:,1),xGmm(:,2),[],'filled');
                %scatter3(xGmm(:,1),xGmm(:,2),compIdx(:,i+1),[],compIdx(:,i+1),'filled');
                %colormap(jet(sDistParams.GMModel.NumComponents));
                %colorbar('TickLabelInterpreter', 'latex');
                grid on;
                xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
                ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
                view(2);
                set(gca,'FontSize', 14);
                xlim([xMin(1), xMax(1)])
                ylim([xMin(2), xMax(2)])
            elseif dim == 3
                scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[], 'filled');
                %scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[],compIdx(:,i+1), 'filled');
                %colormap(jet(sDistParams.GMModel.NumComponents));
                %colorbar('TickLabelInterpreter', 'latex');
                xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
                ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
                zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
                view(30,75);
                set(gca,'FontSize', 14);
                xlim([xMin(1), xMax(1)])
                ylim([xMin(2), xMax(2)])
                zlim([xMin(3), xMax(3)])
            end
            title(plt2Title, 'Interpreter', 'latex', 'FontSize', 14)
        end
        vAx(1) = nexttile(1);
        UpdateCursorDataTip(fig, vAx, compIdx);
        
    else
        xMax = max(x);
        xMin = min(x);
    end

    %% Dataset
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
        xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
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
    title(strcat(pltTitle, " (", actualDataDist, ")"), 'Interpreter', 'latex', 'FontSize', 14)
    
    %% Size
    if strcmp(windowStyle, 'normal')
        x0     = 400;
        y0     = 400;
        height = 400;
        width  = 600 + 600*(exist('b_plotDistModel', 'var') && b_plotDistModel);
        set(gcf,'Position', [x0 y0 width height])
    end
    set(0,'DefaultFigureWindowStyle',prevWindowStyle)
end
%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder') && exist('fig','var')
    if exist('b_plotDistModel', 'var') && b_plotDistModel
        figName = 'dataset_vs_dist';
    else
        figName = 'dataset';
    end
    set(fig,'renderer','Painters')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end