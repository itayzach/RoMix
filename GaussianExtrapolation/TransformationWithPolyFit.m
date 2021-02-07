%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = false;
%% Random data
verticesPDF = 'SwissRoll'; % 'Gaussian' / 'Uniform' / 'TwoMoons' / 'TwoSpirals'
n = 1000;
dim = 3;
if strcmp(verticesPDF, 'Uniform')
    xMax = 1*ones(1,dim);
    xMin = -1*ones(1,dim);
%     n = round(sqrt(n))^2;
%     [X1, X2] = meshgrid(linspace(xMin(1),xMax(1),sqrt(n)),linspace(xMin(2),xMax(2),sqrt(n)));
%     xTrain = [X1(:) X2(:)];
    xTrain = (xMax - xMin).*rand(n,dim) + xMin;
    omega = 0.1;
elseif strcmp(verticesPDF, 'Gaussian')
    sigma = 1;
    mu = 0;
    xMax = 3*sigma*ones(1,dim);
    xMin = -3*sigma*ones(1,dim);
    xTrain = sigma*randn(n,dim) + mu;
    omega = 1;
elseif strcmp(verticesPDF, 'TwoMoons')
    b_loadTwoMoonsMatFile = false;
    N = 5000;
    sDataset = GenerateTwoMoonsDataset(n, N, b_loadTwoMoonsMatFile);
    xTrain = sDataset.x;
    n = length(xTrain);
    N = length(xTest);
    xMax = max(xTrain);
    xMin = min(xTrain);
    omega = 0.3;
    assert(dim == 2);
elseif strcmp(verticesPDF, 'TwoSpirals')
    N = 5000;
    sTwoSpirals = GenerateTwoSpiralsDataset(n, N);
    xTrain = sTwoSpirals.x;
    xMax = max(xTrain);
    xMin = min(xTrain);
    omega = 0.3;
    assert(dim == 2);
elseif strcmp(verticesPDF, 'SwissRoll')
    xTrain = GenerateSwissRoll(n);
    xMax = max(xTrain);
    xMin = min(xTrain);
    omega = 0.3;
    assert(dim == 3);
    N = 5000;
    sDataset.xt = GenerateSwissRoll(N);
else
    error('invalid verticesPDF');
end

if dim <= 2
    PlotHistogram(sSimParams, xTrain, verticesPDF, 'Histogram of X', false);
end
%% Generate graph
dist = pdist2(xTrain, xTrain);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
M = 20;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:4;

if dim <= 3
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, V, lambda, 'Numeric', [], [], 'Eigenvectors of ${\bf W}$')
end
%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
pCdfDegree = 10;
invpCdfDegree = 10;
[xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = PolyfitEstCdf(xTrain, n*ones(1,dim), pCdfDegree, invpCdfDegree);

%% Transform to Gtilde
muTilde = 0*ones(1,dim);
sigmaTilde = 1*eye(dim);
xTildeTrain = T(pCdf, true, muTilde, sigmaTilde, xTrain);

if dim == 2
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
end
PlotHistogram(sSimParams, xTildeTrain, 'Gaussian', 'Histogram of X', false);
%% Demonstrate T
for d = 1:dim
    PlotPolyCdfDemonstration1(xMin(d), xMax(d), pCdf(d,:), xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), muTilde(d), sigmaTilde(d,d));
    PlotPolyCdfDemonstration2(xMin(d), xMax(d), pCdf(d,:), invpCdf(d,:), muTilde(d), sigmaTilde(d,d));
end

%% Build G tilde
omegaTilde = 0.3;
distTilde = pdist2(xTildeTrain, xTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));
%% Calculate analytic eigenfunctions of W_tilde, and numeric eigenvectors for comparison
MTilde = 20;
[PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

if dim <= 3
    figTitle = 'Eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1)-1, vInd(end)-1, VTilde, lambdaNumericTilde, 'Numeric', [], [], figTitle)
    figTitle = 'Eigenfunctions of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1)-1, vInd(end)-1, PhiTilde, lambdaAnalyticTilde, 'Analytic', [], [], figTitle)
end
%% Functional maps
C = pinv(PhiTilde)*V;

figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

%% V in terms of Phi_tilde
VRec = PhiTilde*C;

if dim <= 3
    figTitle = ['${\bf V}$ in terms of ${\bf \tilde{\Phi}} \quad ({\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C})$'];
    figName = 'VRec';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, VRec, lambda, 'Numeric', [], [], figTitle, figName)
end
%% Interpolate
if ~exist('N', 'var')
    N = 5000;
    if dim == 2
        [X1, X2] = meshgrid(linspace(xMin(1),xMax(1),sqrt(N)),linspace(xMin(2),xMax(2),sqrt(N)));
        xInt = [X1(:) X2(:)];
    elseif dim == 1
        xInt = linspace(xMin, xMax, N)';
    end
else
    xInt = sDataset.xt;
end

xTildeInt = T(pCdf, true, muTilde, sigmaTilde, xInt);
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);

if dim <= 3
    if dim == 3
        figTitle = ['Eigenfunctions of $\tilde{{\bf W}}$ on the entire plane' newline ...
            '$\tilde{\omega}$ = ' num2str(omegaTilde, '%.2f') ...
            '; $\tilde{\Sigma} =  \left[ {\matrix{ ',num2str(sigmaTilde(1,1)),' & ',  num2str(sigmaTilde(1,2)),' & ',  num2str(sigmaTilde(1,3)),...
            ' \cr ', num2str(sigmaTilde(2,1)) , ' &  ',  num2str(sigmaTilde(2,2)),' & ',  num2str(sigmaTilde(2,3)),...
            ' \cr ', num2str(sigmaTilde(3,1)) , ' &  ',  num2str(sigmaTilde(3,2)),' & ',  num2str(sigmaTilde(3,3)),' } }  \right]$', ...
            '; $\tilde{\mu} =  \left[ {\matrix{ ',num2str(muTilde(1)), ...
            ' \cr ', num2str(muTilde(3)), ...
            ' \cr ', num2str(muTilde(2)),' } }  \right]$' ];
    elseif dim == 2
        figTitle = ['Eigenfunctions of $\tilde{{\bf W}}$ on the entire plane' newline ...
            '$\tilde{\omega}$ = ' num2str(omegaTilde, '%.2f') ...
            '; $\tilde{\Sigma} =  \left[ {\matrix{ ',num2str(sigmaTilde(1,1)),' & ',  num2str(sigmaTilde(1,2)),...
            ' \cr ', num2str(sigmaTilde(2,1)) , ' &  ',  num2str(sigmaTilde(2,2)),' } }  \right]$', ...
            '; $\tilde{\mu} =  \left[ {\matrix{ ',num2str(muTilde(1)), ...
            ' \cr ', num2str(muTilde(2)),' } }  \right]$' ];
    elseif dim == 1
        figTitle = ['Eigenfunctions of $\tilde{{\bf W}}$ on the entire axis' newline ...
        '$\tilde{\omega}$ = ' num2str(omegaTilde, '%.2f') ...
        '; $\tilde{\sigma}$ = ' num2str(sigmaTilde, '%.2f') ...
        '; $\tilde{\mu}$ = ' num2str(muTilde, '%.2f') ];
    end
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeInt, [], vInd(1)-1, vInd(end)-1, PhiTildeInt, lambdaAnalyticTilde, 'Numeric', [], [], figTitle)
end

% xIntInvT = invT(invpCdf, muTilde, sigmaTilde, xTildeInt);
VInt = PhiTildeInt*C;

interpRatio = N/n;
VIntRenormed = sqrt(interpRatio)*VInt;

if dim <= 3
    figTitle = ['Interpolated eigenvectors of ${\bf W}$'];% newline '${\bf V}_{{\bf int}} = \sqrt{\frac{N}{n}}{\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'];
    figName = 'VInt';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1)-1, vInd(end)-1, VIntRenormed, lambda, 'Numeric', [], [], figTitle, figName)
end

