function PlotAicBic(sPlotParams, sPreset, sDataset, sDistParams)

vGmmNumCompAicBic = sPreset.vGmmNumCompAicBic;
gmmMaxIter = sPreset.gmmMaxIter;
x = sDataset.sData.x;
gmmRegVal = sPreset.gmmRegVal;

assert(ismember(sDistParams.GMModel.NumComponents, vGmmNumCompAicBic))
vGMModels = cell(size(vGmmNumCompAicBic(:),1),1);
vAIC = zeros(size(vGmmNumCompAicBic(:),1),1);
vBIC = zeros(size(vGmmNumCompAicBic(:),1),1);
i = 1;
for nComp = vGmmNumCompAicBic
    if nComp == sDistParams.GMModel.NumComponents
        vGMModels{i} = sDistParams.GMModel;
        fprintf('Taking fitgmdist from sDistParams with %d components... ', nComp)
    else
        fprintf('Running fitgmdist with %d components... ', nComp)
        converged = false;
        numAttempts = 0;
        warning('off','stats:gmdistribution:FailedToConverge');
        totalAttempts = 10;
        while ~converged && numAttempts < totalAttempts
            options = statset('MaxIter',gmmMaxIter);
            vGMModels{i} = fitgmdist(x, nComp, 'RegularizationValue', gmmRegVal, 'Options', options);
            converged = vGMModels{i}.Converged;
            numAttempts = numAttempts + 1;
            fprintf('%d... ', numAttempts)
        end
        assert(converged, 'GMM couldn''t converge...')
        fprintf('Done after %d iterations with AIC = %.2f, BIC = %.2f\n', vGMModels{i}.NumIterations, vGMModels{i}.AIC, vGMModels{i}.BIC)
    end
    vAIC(i) = vGMModels{i}.AIC;
    vBIC(i) = vGMModels{i}.BIC;
    i = i + 1;
end

%% Plot
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')
fig = figure; 
%plot(vGmmNumCompAicBic, vAIC, '-o', 'DisplayName', 'AIC', 'LineWidth',2); 
%hold on; 
plot(vGmmNumCompAicBic, vBIC, '-o', 'DisplayName', 'BIC', 'LineWidth',2);
%legend();
xlim([min(vGmmNumCompAicBic), max(vGmmNumCompAicBic)])
ylabel('BIC', 'Interpreter','latex','FontSize',14)
xlabel('$k$', 'Interpreter','latex','FontSize',14)
set(gca,'FontSize', 14);
set(gca,'YTickLabel',[])
%set(gca, 'YScale', 'log');
%set(gca, 'YTickLabel', {'$10^{0}$', '$10^{1}$'})
%set(gca, 'YTick', [1 10])
%% Size
x0     = 10;
y0     = 50;
height = 200;
width  = 600;
set(gcf,'Position', [x0 y0 width height])

%% Save
SaveFigure(sPlotParams, fig, 'BIC', {'epsc', 'png'});
set(0,'DefaultFigureWindowStyle',windowStyle)
end