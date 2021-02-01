%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = false;
%% Random data
verticesPDF = 'Uniform'; % 'Gaussian' / 'Uniform'
n = 2000;
dim = 2;
if strcmp(verticesPDF, 'Uniform')
    xMax = [1 1];
    xMin = [-1 -1];
%     n = round(sqrt(n))^2;
%     [X1, X2] = meshgrid(linspace(xMin(1),xMax(1),sqrt(n)),linspace(xMin(2),xMax(2),sqrt(n)));
%     xTrain = [X1(:) X2(:)];
    xTrain = (xMax - xMin).*rand(n,dim) + xMin;
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
omega = 0.3;
dist = pdist2(xTrain, xTrain);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
M = 50;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:10;

PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, V, lambda, 'Numeric', [], [], 'Eigenvectors of ${\bf W}$')
%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
assert(dim == 2, 'following function only works for 2D')
[estCdf_xTrainGrid, xTrainGrid, estMarginalCdf_xTrain] = ecdf2(xTrain, [n n]);

pCdfDegree = 2;
invpCdfDegree = 2;
pCdf{1} = polyfit(xTrainGrid{1}, estMarginalCdf_xTrain{1}, pCdfDegree); % to have analytic expression for the cdf
invpCdf{1} = polyfit(estMarginalCdf_xTrain{1}, xTrainGrid{1}, invpCdfDegree); % to have analytic expression for the cdf
pCdf{2} = polyfit(xTrainGrid{2}, estMarginalCdf_xTrain{2}, pCdfDegree); % to have analytic expression for the cdf
invpCdf{2} = polyfit(estMarginalCdf_xTrain{2}, xTrainGrid{2}, invpCdfDegree); % to have analytic expression for the cdf


%% Transform to Gtilde
muTilde = [0 0];
sigmaTilde = [1 0; 0 1];
xTildeTrain(:,1) = T(pCdf{1}, true, muTilde(1), sigmaTilde(1,1), xTrain(:,1));
xTildeTrain(:,2) = T(pCdf{2}, true, muTilde(2), sigmaTilde(2,2), xTrain(:,2));

figure('Name', 'x vs. xTilde');
subplot(1,2,1)
scatter(xTrain(:,1),xTrain(:,2),50,1:n,'filled');
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
subplot(1,2,2)
scatter(xTildeTrain(:,1),xTildeTrain(:,2),50,1:n,'filled');
xlabel('$\tilde{x}_1$', 'interpreter', 'latex', 'FontSize', 16);
ylabel('$\tilde{x}_2$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

PlotHistogram(sSimParams, xTildeTrain, 'Gaussian', 'Histogram of X', false);
%% Demonstrate T (1/2)
PlotPolyCdfDemonstration1(xMin(1), xMax(1), pCdf{1}, xTrainGrid{1}, estMarginalCdf_xTrain{1}, muTilde(1), sigmaTilde(1,1));
PlotPolyCdfDemonstration1(xMin(2), xMax(2), pCdf{2}, xTrainGrid{2}, estMarginalCdf_xTrain{2}, muTilde(2), sigmaTilde(2,2));

%% Demonstrate T (2/2)
PlotPolyCdfDemonstration2(xMin(1), xMax(1), pCdf{1}, invpCdf{1}, muTilde(1), sigmaTilde(1,1));
PlotPolyCdfDemonstration2(xMin(2), xMax(2), pCdf{2}, invpCdf{2}, muTilde(2), sigmaTilde(2,2));

%% Build G tilde
omegaTilde = 0.3;
distTilde = pdist2(xTildeTrain, xTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));
%% Calculate analytic eigenfunctions of W_tilde, and numeric eigenvectors for comparison
MTilde = 50;
[PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

figTitle = 'Eigenvectors $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1)-1, vInd(end)-1, VTilde, lambdaNumericTilde, 'Numeric', [], [], figTitle)
figTitle = 'Eigenfunctions $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1)-1, vInd(end)-1, PhiTilde, lambdaAnalyticTilde, 'Analytic', [], [], figTitle)
%% Functional maps
C = pinv(PhiTilde)*V;
pinvC = pinv(C);

figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

figure('Name', 'pinv(C)');
imagesc(pinvC);
colorbar();
title('${\bf C}^\dagger$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

%% V in terms of Phi_tilde
VRec = PhiTilde*C;

figTitle = ['${\bf V}$ in terms of ${\bf \tilde{\Phi}}$' newline '${\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C}$'];
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, VRec, lambda, 'Numeric', [], [], figTitle)
%% Interpolate
N = 5000;
[X1, X2] = meshgrid(linspace(xMin(1),xMax(1),sqrt(N)),linspace(xMin(2),xMax(2),sqrt(N)));
xInt = [X1(:) X2(:)];
xTildeInt(:,1) = T(pCdf{1}, true, muTilde(1), sigmaTilde(1,1), xInt(:,1));
xTildeInt(:,2) = T(pCdf{2}, true, muTilde(2), sigmaTilde(2,2), xInt(:,2));
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);

figTitle = ['Eigenfunctions of $\tilde{{\bf W}}$ on the entire axis' ];%newline ...
%     '$\tilde{\omega}$ = ' num2str(omegaTilde, '%.2f') ...
%     '; $\tilde{\sigma}$ = ' num2str(sigmaTilde, '%.2f') ...
%     '; $\tilde{\mu}$ = ' num2str(muTilde, '%.2f') ]; 
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeInt, [], vInd(1)-1, vInd(end)-1, PhiTildeInt, lambdaAnalyticTilde, 'Numeric', [], [], figTitle)


xIntInvT(:,1) = invT(invpCdf{1}, muTilde(1), sigmaTilde(1,1), xTildeInt(:,1));
xIntInvT(:,2) = invT(invpCdf{2}, muTilde(2), sigmaTilde(2,2), xTildeInt(:,2));
VInt = PhiTildeInt*C;

interpRatio = N/n;
VIntRenormed = sqrt(interpRatio)*VInt;

figTitle = ['Interpolated eigenvectors of ${\bf W}$' newline '${\bf V}_{{\bf int}} = \sqrt{\frac{N}{n}}{\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'];
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xIntInvT, [], vInd(1)-1, vInd(end)-1, VIntRenormed, lambda, 'Numeric', [], [], figTitle)


%% ecdf2d
function [estCdf_xGrid, xGrid, estMarginalCdf] = ecdf2(x, bins)
[hist_xGrid, xGrid] = hist3(x, 'Nbins', bins);
n = length(x);
pdf_xGrid = (1/n)*hist_xGrid;
estCdf_xGrid = cumsum(cumsum(pdf_xGrid,1),2);  % CDF_x1x2
estMarginalCdf{1} = cumsum(sum(pdf_xGrid,2));  % CDF_x1
estMarginalCdf{2} = cumsum(sum(pdf_xGrid,1))'; % CDF_x2
end