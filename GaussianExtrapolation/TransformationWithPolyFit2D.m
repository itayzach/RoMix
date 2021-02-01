%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = true;
%% Random data
verticesPDF = 'Uniform'; % 'Gaussian' / 'Uniform'
n = 1000;
dim = 2;
if strcmp(verticesPDF, 'Uniform')
    xMax = 10;
    xMin = -10;
%     n = round(sqrt(n))^2;
%     [X1, X2] = meshgrid(linspace(xMin,xMax,sqrt(n)));
%     xTrain = [X1(:) X2(:)];
    xTrain = (xMax - xMin)*rand(n,dim) + xMin;
elseif strcmp(verticesPDF, 'Gaussian')
    sigma = 10;
    mu = 0;
    xMax = 3*sigma;
    xMin = -3*sigma;
    xTrain = sigma*randn(n,dim) + mu;
else
    error('invalid verticesPDF');
end

PlotHistogram(sSimParams, xTrain, verticesPDF, 'Histogram of X', false);
%% Generate graph
omega = 3;
dist = pdist2(xTrain, xTrain);
W = exp(-dist.^2/(2*omega^2));
G = gsp_graph(W, xTrain);

%% Calculate (numeric) eigenvectors of G
M = 20;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:5;

sSimParams.PlotEigenFuncsM = M;
% PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, V, lambda, 'Numeric', [], G, 'Eigenvectors of ${\bf W}$')
%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
assert(dim == 2, 'following function only works for 2D')
[estCdf_xTrainGrid, xTrainGrid, estCdf_xTrain1, estCdf_xTrain2] = ecdf2(xTrain, [n n]);

figure;
surf(xTrainGrid{1}, xTrainGrid{2}, estCdf_xTrainGrid', 'EdgeColor', 'none');
% hold on;
figure;
plot(xTrainGrid{1}', estCdf_xTrain1, '.');
hold on;
plot(xTrainGrid{2}', estCdf_xTrain2, '.');


pCdfDegree = 10;
[X_c, Y_c] = meshgrid(xTrainGrid{1}, xTrainGrid{2});
xTrainMeshGrid = [X_c(:), Y_c(:)];
sPolyCdf = polyfitn(xTrainMeshGrid,estCdf_xTrainGrid(:),pCdfDegree); % to have analytic expression for the cdf

estCdf_xTrain = polyvaln(sPolyCdf, xTrain);
estCdf_xTrain(estCdf_xTrain<0) = eps; % saturate

figure; 
surf(xTrainGrid{1}, xTrainGrid{2}, estCdf_xTrainGrid', 'EdgeColor', 'none');
hold on;
plot3(xTrain(:,1), xTrain(:,2), estCdf_xTrain, 'b.');
legend('${\bf eCDF}_{{\bf X}}(x_{{\bf train}})$', '$\hat{F}_{X}(x_{{\bf test}})$', ...
        'location', 'northeast', 'interpreter', 'latex', 'FontSize', 14)
title(['Est. vs. poly CDF (degree ' num2str(pCdfDegree) ')'], 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
view(80,25);

% figure; 
% This is just a reference
% [X1, X2] = meshgrid(-4:0.05:4,-4:0.05:4);
% xGrid = [X1(:) X2(:)];
% mvncdfval = mvncdf(xGrid, [0 0], [1 0; 0 1]);
% surf(-4:0.05:4,-4:0.05:4,reshape(mvncdfval, size(X1,1), size(X1,2)), 'EdgeColor', 'none')
% title('CDF of Multivariate Gaussian', 'interpreter', 'latex', 'FontSize', 16);
% set(gca,'FontSize', 14);
% view(80,25);

muTilde = 0;
sigmaTilde = 1;
xTilde1 = icdf('Normal',estCdf_xTrain1,muTilde,sigmaTilde);
xTilde2 = icdf('Normal',estCdf_xTrain2,muTilde,sigmaTilde);
xTilde = [xTilde1 xTilde2];

% figure;
% plot3(xTilde1,xTilde2,estCdf_xTrain,'.')

figure; 
subplot(1,2,1)
histfit(xTilde1(xTilde1 < Inf),100);
subplot(1,2,2)
histfit(xTilde2(xTilde2 < Inf),100);
%% Transform to Gtilde
muTilde = 0*ones(dim,1);
sigmaTilde = 1*eye(dim); % independent
xTildeTrain = Td(sPolyCdf, true, muTilde, sigmaTilde, xTrain);

%% Demonstrate T (1/2)
% xTestGrid = linspace(xMin,xMax,2000)';
% polyCdf_xTestGrid = polyval(pCdf, xTestGrid);
% b_saturate = false;
% if b_saturate && (any(polyCdf_xTestGrid > 1) || any(polyCdf_xTestGrid < 0))
%     warning('CDF must be in [0,1], saturating...');
%     polyCdf_xTestGrid(polyCdf_xTestGrid > 1) = 1-eps; % saturate
%     polyCdf_xTestGrid(polyCdf_xTestGrid < 0) = eps; % saturate
% end
% 
% xTildeTestGrid = icdf('Normal',polyCdf_xTestGrid,muTilde,sigmaTilde);
% 
% % Generate the polynomial title
% pCdfStr = [];
% for p = 1:pCdfDegree+1
%     if abs(pCdf(p)) < 1e-4
%         continue;
%     end
%     if p > 1 && pCdf(p) > 0 && ~isempty(pCdfStr)
%         pCdfStr = strcat(pCdfStr,'+');
%     end
%     if p == pCdfDegree+1
%         pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree+1),'%.5f'));
%     elseif p == pCdfDegree
%         pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree),'%.5f'),'x');
%     else
%         pCdfStr = strcat(pCdfStr, num2str(pCdf(p),'%.5f'),'x^{',num2str(pCdfDegree-p+1),'}');
%     end
%     
% end
% figure('Name', 'Demonstrate T (1/2)');
% subplot(2,2,1)
%     plot(xTrainGrid, estCdf_xTrainGrid,'.');
%     hold on;
%     plot(xTestGrid, polyCdf_xTestGrid, '.');
%     ylim([0 1])
%     title(['Est. vs. poly CDF (degree ' num2str(pCdfDegree) ')'], 'interpreter', 'latex', 'FontSize', 16);
%     xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
%     ylabel('$\hat{F}_{X}(x)$', 'interpreter', 'latex', 'FontSize', 16);
%     legend('${\bf eCDF}_{{\bf X}}(x_{{\bf train}})$', '$\hat{F}_{X}(x_{{\bf test}})$', ...
%         'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
%     set(gca,'FontSize', 14);
% subplot(2,2,2)
%     plot(xTildeTestGrid,polyCdf_xTestGrid,'.')
%     title(['$\tilde{x} = T(x) = F_{\tilde{X}}^{-1}(\hat{F}_{X}(x))$' newline '(flipped axes)'], 'interpreter', 'latex', 'FontSize', 16);
%     ylabel('$\hat{F}_{X}(x)$', 'interpreter', 'latex', 'FontSize', 16);
%     xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
%     set(gca,'FontSize', 14);
% subplot(2,2,3)
%     histogram(xTestGrid ,100);
%     title('Histogram of $x_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
%     set(gca,'FontSize', 14);
% subplot(2,2,4)
%     histfit(xTildeTestGrid ,100);
%     title('Histogram of $\tilde{x}_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
%     set(gca,'FontSize', 14);
% sgtitle(['$\hat{F}_{X}(x) = ' pCdfStr '$'], 'interpreter', 'latex', 'FontSize', 16);
% saveas(fig,strcat(outputFolder, filesep, 'fig3_polyfit'), figSaveType);

%% Demonstrate T (2/2)
% nTestPoints = 50;
% xSmallTest = linspace(xMin+0.5,xMax-0.5,nTestPoints)';
% 
% xTestTilde = T(pCdf, true, muTilde, sigmaTilde, xSmallTest);
% xSmallTestEst = invT(invpCdf, muTilde, sigmaTilde, xTestTilde);
% 
% cmap = xSmallTest;
% figure('Name', 'Demonstrate T (2/2)');
% subplot(2,1,1)
%     scatter(xSmallTest, zeros(1,nTestPoints), 100, cmap, 'o')
%     hold on;
%     scatter(xSmallTestEst, zeros(1,nTestPoints), 50, cmap, 'filled')
%     colormap('jet')
%     xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
%     legend('$x_{{\bf test}}$', '$T^{-1}(T(x_{{\bf test}}))$','interpreter', 'latex', 'FontSize', 14);
%     title('Original nodes','interpreter', 'latex', 'FontSize', 16);
%     set(gca,'YTick',[],'FontSize', 14);
% subplot(2,1,2);
%     scatter(xTestTilde, zeros(1,nTestPoints), 50, cmap, 'filled')
%     colormap('jet')
%     xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
%     title('Transformed nodes','interpreter', 'latex', 'FontSize', 16);
%     set(gca,'YTick',[],'FontSize', 14);
% saveas(fig,strcat(outputFolder, filesep, 'fig4_x_to_xTilde'), figSaveType);

%% Build G tilde
omegaTilde = 0.3;
distTilde = pdist2(xTildeTrain, xTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));

%% Calculate analytic eigenfunctions of W_tilde, and numeric eigenvectors for comparison
MTilde = 20;
[PhiTilde, ~] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

figure('Name', 'evecs vs. efuncs of WTilde');
plot(xTildeTrain, VTilde(:,vInd),'o');
hold on
plot(xTildeTrain, PhiTilde(:,vInd),'.');
legend([strcat('$\tilde{v}_{',string(vInd),'}$') strcat('$\tilde{\phi}_',string(vInd),'$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('Eigenfunctions vs. eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Functional maps
C = pinv(PhiTilde)*V;
pinvC = pinv(C);

figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

figure('Name', 'pinv(C)');
imagesc(pinvC);
colorbar();
title('${\bf C}^\dagger$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% V in terms of Phi_tilde
VRec = PhiTilde*C;

% figure('Name', 'V in terms of PhiTilde');
% plot(xTrain, V(:,vInd),'o');
% hold on
% plot(xTrain, VRec(:,vInd),'.');
% legend([strcat('$v_{',string(vInd),'}$') strcat('$v_{{\bf rec},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
% title(['${\bf V}$ in terms of ${\bf \tilde{\Phi}}$' newline '${\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); 
% set(gca,'FontSize', 14);


PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], 0, M-1, VRec, lambda, 'Numeric', [], G, 'Eigenvectors of ${\bf W}$')
%% Interpolate
N = 5000;
xInt = linspace(xMin+0.5, xMax-0.5, N)';
xTildeInt = Td(pCdf, true, muTilde, sigmaTilde, xInt);
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);

figure('Name', 'Eigenfunctions of WTilde on the entire axis');
plot(xTildeInt, PhiTildeInt(:,vInd),'LineWidth',2);
legend(strcat('$\tilde{\phi}_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Eigenfunctions of $\tilde{{\bf W}}$ on the entire axis' newline ...
    '$\tilde{\omega}$ = ' num2str(omegaTilde, '%.2f') ...
    '; $\tilde{\sigma}$ = ' num2str(sigmaTilde, '%.2f') ...
    '; $\tilde{\mu}$ = ' num2str(muTilde, '%.2f') ], 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

xIntInvT = invT(invpCdf, muTilde, sigmaTilde, xTildeInt);
VInt = PhiTildeInt*C;

interpRatio = N/n;
VIntRenormed = sqrt(interpRatio)*VInt;

figure('Name', 'Interpolated evecs of W');
plot(xTrain, V(:,vInd),'o');
hold on
plot(xIntInvT, VIntRenormed(:,vInd),'.');
legend([strcat('$v_{',string(vInd),'}$') strcat('$v_{{\bf int},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Interpolated eigenvectors of ${\bf W}$' newline '${\bf V}_{{\bf int}} = \sqrt{\frac{N}{n}}{\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Td
function xTilde = Td(sPolyCdf, b_saturate, mu, sigma, x)
polyCdf = polyvaln(sPolyCdf, x);
if (any(polyCdf > 1) || any(polyCdf < 0))
    warning('CDF must be in [0,1]...');
    if b_saturate
        polyCdf(polyCdf > 0.99) = 0.99; % saturate
        polyCdf(polyCdf < 0.01) = 0.01; % saturate
    end
end
% xTilde = icdf('Normal',polyCdf,mu,sigma);
xTilde = norminv(polyCdf,mu,sigma);
assert(~any(isnan(xTilde)),['xTilde contain NaNs since because of x = ', num2str(x(isnan(xTilde))')]);
end

%% ecdf2d
function [estCdf_xGrid, xGrid, estCdf_x1, estCdf_x2] = ecdf2(x, bins)
[hist_xGrid, xGrid] = hist3(x, 'Nbins', bins);
n = length(x);
% xMax = max(x);
% xMin = min(x);
% dx1 = (xMax(1) - xMin(1))/bins(1);
% dx2 = (xMax(2) - xMin(2))/bins(2);
pdf_xGrid = (1/n)*hist_xGrid;
estCdf_xGrid = cumsum(cumsum(pdf_xGrid,1),2);
figure;
surf(xGrid{1}, xGrid{2}, pdf_xGrid', 'EdgeColor', 'none');
% dx1 = [0 diff(sort(xGrid{1}))]';
% dx2 = [0 diff(sort(xGrid{2}))]';
estCdf_x1 = cumsum(sum(pdf_xGrid,2));
estCdf_x2 = cumsum(sum(pdf_xGrid,1))';
end