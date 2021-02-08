function [] = PlotPolyCdfDemonstration1(xMin, xMax, pCdf, xTrainGrid, estCdf_xTrainGrid, muTilde, sigmaTilde)
pCdfDegree = length(pCdf)-1;
xTestGrid = linspace(xMin,xMax,2000)';
polyCdf_xTestGrid = polyval(pCdf, xTestGrid);
b_saturate = false;
if b_saturate && (any(polyCdf_xTestGrid > 1) || any(polyCdf_xTestGrid < 0))
    warning('CDF must be in [0,1], saturating...');
    polyCdf_xTestGrid(polyCdf_xTestGrid > 1) = 1-eps; % saturate
    polyCdf_xTestGrid(polyCdf_xTestGrid < 0) = eps; % saturate
end

xTildeTestGrid = icdf('Normal',polyCdf_xTestGrid,muTilde,sigmaTilde);

% Generate the polynomial title
pCdfStr = '\hat{F}_{X}(x) = ';
pCdfCells = {};
for p = 1:pCdfDegree+1
    if abs(pCdf(p)) < 1e-4
        continue;
    end
    if pCdfDegree > 5
        b_splitToTwoLines = true;
    else
        b_splitToTwoLines = false;
    end
    if p > 1 && pCdf(p) > 0 && ~isempty(pCdfStr)
        pCdfStr = strcat(pCdfStr,'+');
    end
    if p == pCdfDegree+1
        pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree+1),'%.5f'));
    elseif p == pCdfDegree
        pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree),'%.5f'),'x');
    else
        pCdfStr = strcat(pCdfStr, num2str(pCdf(p),'%.5f'),'x^{',num2str(pCdfDegree-p+1),'}');
        if p == floor(pCdfDegree/2) && b_splitToTwoLines
            pCdfCells{end+1} = [ '$' pCdfStr '$'];
            pCdfStr = [];
        end
    end
end
pCdfCells{end+1} = [ '$' pCdfStr '$'];

fig = figure('Name', 'Demonstrate T #1');
subplot(2,2,1)
    plot(xTrainGrid, estCdf_xTrainGrid,'.');
    hold on;
    plot(xTestGrid, polyCdf_xTestGrid, '.');
    ylim([0 1])
    title(['Est. vs. poly CDF (degree ' num2str(pCdfDegree) ')'], 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$\hat{F}_{X}(x)$', 'interpreter', 'latex', 'FontSize', 16);
    legend('${\bf eCDF}_{{\bf X}}(x_{{\bf train}})$', '$\hat{F}_{X}(x_{{\bf test}})$', ...
        'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(2,2,2)
    plot(xTildeTestGrid,polyCdf_xTestGrid,'.')
    title(['$\tilde{x} = T(x) = F_{\tilde{X}}^{-1}(\hat{F}_{X}(x))$' newline '(flipped axes)'], 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$\hat{F}_{X}(x)$', 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(2,2,3)
    histogram(xTestGrid ,100);
    title('Histogram of $x_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(2,2,4)
    histfit(xTildeTestGrid ,100);
    title('Histogram of $\tilde{x}_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
sgtitle(pCdfCells ,'interpreter', 'latex', 'FontSize', 16);
% saveas(fig,strcat(outputFolder, filesep, 'fig3_polyfit'), figSaveType);
end

