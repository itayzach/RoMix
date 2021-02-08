%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = false;
%% Random data
dim = 3;
nComponents = 1;
n = 800;
N = 5000;
verticesPDF = 'SwissRoll'; % 'Gaussian' / 'Uniform' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N);
dim = sDataset.dim;
xTrain = sDataset.sData.x;
xMax = sDataset.xMax;
xMin = sDataset.xMin;

%% Histogram
if dim <= 2
    PlotHistogram(sSimParams, xTrain, verticesPDF, 'Histogram of X', false);
end
%% Generate graph
omega = sDataset.recommendedOmega;
dist = pdist2(xTrain, xTrain);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
M = 50;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:4;

figTitle = 'Eigenvectors of ${\bf W}$';
figName = 'V';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1)-1, vInd(end)-1, V, lambda, [], [], figTitle, figName, 'v')
%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
pCdfDegree = 10;
invpCdfDegree = 10;
[xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = PolyfitEstCdf(xTrain, n*ones(1,dim), pCdfDegree, invpCdfDegree);

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

figTitle = ['${\bf V}$ in terms of ${\bf \tilde{\Phi}} \quad ({\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C})$'];
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

figTitle = ['Interpolated eigenvectors of ${\bf W}$'];% newline '${\bf V}_{{\bf int}} = \sqrt{\frac{N}{n}}{\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'];
figName = 'VInt';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1)-1, vInd(end)-1, VIntRenormed, lambda, [], [], figTitle, figName, 'v^{{\bf int}}')


function figTitle = GetInterpFigTitle(dim, omegaTilde, sigmaTilde, muTilde)
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

end


