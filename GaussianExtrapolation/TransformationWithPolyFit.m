%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = false;
sSimParams.b_showLaplacianEvecs       = false;
sSimParams.b_plotWeights              = false;
%% Random data
dim = 2;
nComponents = 1;
n = 1000;
N = 5000;
verticesPDF = 'SwissRoll'; % 'SwissRoll' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N);
dim = sDataset.dim;
xTrain = sDataset.sData.x;
xMax = sDataset.xMax;
xMin = sDataset.xMin;
n = length(sDataset.sData.x);
N = length(sDataset.sData.xt);

%% Histogram
if dim <= 2
    PlotHistogram(sSimParams, xTrain, verticesPDF, 'Histogram of X', false);
end
PlotDataset(sSimParams, xTrain, verticesPDF, 'Training set');
%% Generate graph
omega = sDataset.recommendedOmega;
dist = pdist2(xTrain, xTrain);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
M = 50;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
if dim > 1
    vInd = 1:12;
else
    vInd = 1:4;
end

figTitle = [ 'Eigenvectors of ${\bf W}$ with $\omega = ' num2str(omega) '\quad n = ' num2str(n) '$'];
figName = 'V';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, ...
    V(:,vInd), lambda, [], [], figTitle, figName, 'v')

%% Plot W on first 4 nodes
if sSimParams.b_plotWeights
    figTitle = [ '${\bf W}(:,i)$ for first $i=1:4$ nodes with $\omega = ' num2str(omega) ...
        '\quad n = ' num2str(n) '$'];
    figName = 'W';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], 0, 3, W(:,1:4), lambda, [], [], ...
        figTitle, figName, 'w')
end
%% Laplacian eigenvectors
if sSimParams.b_showLaplacianEvecs
    D = sum(W,1);
    L = D - W;
    [U, LambdaL] = eigs(L,M);
    lambdaL = diag(LambdaL);
    U = real(U);
    U = FlipSign(V, U);
    figTitle = 'Eigenvectors of ${\bf L}$';
    figName = 'U';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, ...
        U(:,vInd), lambdaL, [], [], figTitle, figName, 'v')
end
%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
pCdfDegree = 10;
invpCdfDegree = 10;
bins = min(700, n);
[xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = PolyfitEstCdf(xTrain, bins, pCdfDegree, invpCdfDegree);

%% Transform to Gtilde
muTilde = mean(xTrain);
sigmaTilde = diag(std(xTrain));
xTildeTrain = T(pCdf, true, muTilde, sigmaTilde, xTrain);

if dim <= 2
    PlotHistogram(sSimParams, xTildeTrain, 'Gaussian', 'Histogram of $\tilde{X}$', false);
end
%% Demonstrate T
for d = 1:dim
    PlotPolyCdfDemonstration1(xMin(d), xMax(d), pCdf(d,:), xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), muTilde(d), sigmaTilde(d,d));
    PlotPolyCdfDemonstration2(xMin(d), xMax(d), pCdf(d,:), invpCdf(d,:), muTilde(d), sigmaTilde(d,d));
end
if dim == 2
    PlotPolyCdfDemonstration3_2D(xTrain,xTildeTrain)
end
%% Build G tilde and numeric eigenvectors for comparison
omegaTilde = omega;
MTilde = 50;
distTilde = pdist2(xTildeTrain, xTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);
%% Calculate analytic eigenfunctions of W_tilde
[PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);

figTitle = 'Eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
figName = 'VTilde';
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1)-1, vInd(end)-1, VTilde, lambdaNumericTilde, [], [], figTitle, figName, '\tilde{v}')
figTitle = 'Eigenfunctions of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
figName = 'PhiTilde';
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1)-1, vInd(end)-1, PhiTilde, lambdaAnalyticTilde, [], [], figTitle, figName, '\tilde{\phi}')

%% Functional maps
C = pinv(PhiTilde)*V;

figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

%% V in terms of Phi_tilde
VRec = PhiTilde*C;

figTitle = '${\bf V}$ in terms of ${\bf \tilde{\Phi}} \quad ({\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C})$';
figName = 'VRec';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, VRec, lambda, [], [], figTitle, figName, 'v^{{\bf rec}}')

b_plotErrVsNodeInd = true;
PlotEigenDiffs(sSimParams, sDataset, [], vInd(1)-1, vInd(end)-1, VRec, V, 'VRec_vs_V', 'v^{{\bf rec}}', 'v', b_plotErrVsNodeInd)
%% Interpolate
if isempty(sDataset.sData.xt)
    N = 5000;
    if dim == 2
        [X1, X2] = meshgrid(linspace(xMin(1),xMax(1),sqrt(N)),linspace(xMin(2),xMax(2),sqrt(N)));
        xInt = [X1(:) X2(:)];
    elseif dim == 1
        xInt = linspace(xMin, xMax, N)';
    end
else
    xInt = sDataset.sData.xt;
end

xTildeInt = T(pCdf, true, muTilde, sigmaTilde, xInt);
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);

figName = 'PhiTildeInt';
figTitle = GetInterpFigTitle(dim, omegaTilde, sigmaTilde, muTilde);
PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeInt, [], vInd(1)-1, vInd(end)-1, PhiTildeInt, lambdaAnalyticTilde, [], [], figTitle, figName, '\tilde{\phi}^{{\bf int}}')

% xIntInvT = invT(invpCdf, muTilde, sigmaTilde, xTildeInt);
VInt = PhiTildeInt*C;

interpRatio = N/n;
VIntRenormed = sqrt(interpRatio)*VInt;

figTitle = 'Interpolated eigenvectors of ${\bf W}$ - Ours';
figName = 'VInt';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1)-1, vInd(end)-1, VIntRenormed, lambda, [], [], figTitle, figName, 'v^{{\bf int}}')

%% Interpolate with Nystrom
distN = pdist2(xTrain, xInt);
WN = exp(-distN.^2/(2*omega^2));

VNys = WN.'*V*diag(1./lambda);

figTitle = 'Interpolated eigenvectors of ${\bf W}$ - Nystrom';
figName = 'VNys';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1)-1, vInd(end)-1, VNys, lambda, [], [], figTitle, figName, 'v^{{\bf nys}}')
