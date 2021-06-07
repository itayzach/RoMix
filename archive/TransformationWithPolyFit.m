%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_plotAllEvecs             = false;
sSimParams.b_GSPBoxPlots              = false;
sSimParams.b_plotLaplacianEvecs       = false;
sSimParams.b_plotWeights              = false;
sSimParams.b_plotVRec                 = false;
sSimParams.b_plotTransDemos           = false;
sSimParams.b_plotOrigVsInterpEvecs    = true;

%% Random data
dim = 1;
nComponents = 1;
n = 1000;
N = 5000;
verticesPDF = 'Grid'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
origGraphAdjacency = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
interpMethod = 'AddPoints'; % 'NewPoints' / 'AddPoints'
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod);
dim = sDataset.dim;
xTrain = sDataset.sData.x;
xMax = sDataset.xMax;
xMin = sDataset.xMin;
n = length(sDataset.sData.x);
N = length(sDataset.sData.xt);

b_kde = false;

%% Histogram
if dim <= 2
    PlotHistogram(sSimParams, xTrain, verticesPDF, 'Histogram of X', false);
end
PlotDataset(sSimParams, xTrain, verticesPDF, 'Training set');
%% Generate graph
omega = sDataset.recommendedOmega;
k = 3;
W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega, k);

%% Calculate (numeric) eigenvectors of G
M = 50;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
% sqrtLambdaV = abs(sqrt(lambda))'.*V;
vInd = 0:4;
% vInd = M-5:M-1;

if sSimParams.b_plotOrigVsInterpEvecs || sSimParams.b_plotAllEvecs
    if strcmp(origGraphAdjacency, 'GaussianKernel')
        figTitle = [ 'Eigenvectors of ${\bf W}$ (Gaussian kernel) with $\omega = ' num2str(omega) '\quad n = ' num2str(n) '$'];
    elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
        figTitle = [ 'Eigenvectors of ${\bf W}$ (k-NN) with $k = ' num2str(k) '\quad n = ' num2str(n) '$'];
    end
    figName = 'V';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1), vInd(end), ...
        V, lambda, [], [], figTitle, figName, 'v')
end
%% Plot W on first 4 nodes
if sSimParams.b_plotWeights
    if strcmp(origGraphAdjacency, 'GaussianKernel')
        figTitle = [ '${\bf W}(:,i)$ (Gaussian kernel) for first $i=1:4$ nodes with $\omega = ' num2str(omega) ...
            '\quad n = ' num2str(n) '$'];
    elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
        figTitle = [ '${\bf W}(:,i)$ (k-NN) for first $i=1:4$ nodes with $k = ' num2str(k) ...
            '\quad n = ' num2str(n) '$'];
    end
    figName = 'W';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], 0, 3, W(:,1:4), lambda, [], [], ...
        figTitle, figName, 'w')
end
%% Laplacian eigenvectors
if sSimParams.b_plotLaplacianEvecs
    D = sum(W,1);
    L = D - W;
    [U, LambdaL] = eigs(L,M);
    lambdaL = diag(LambdaL);
    U = real(U);
    U = FlipSign(V, U);
    figTitle = 'Eigenvectors of ${\bf L}$';
    figName = 'U';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1), vInd(end), ...
        U, lambdaL, [], [], figTitle, figName, 'v')
end
%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
pCdfDegree = 10;
invpCdfDegree = 10;
bins = min(700, n);
[xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = PolyfitEstCdf(xTrain, bins, pCdfDegree, invpCdfDegree, b_kde);

%% Transform to Gtilde
muTilde = mean(xTrain);
sigmaTilde = diag(std(xTrain));
b_saturateT = true;
xTildeTrain = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain);

if dim <= 2
    PlotHistogram(sSimParams, xTildeTrain, 'Gaussian', 'Histogram of $\tilde{X}$', false);
end
%% Demonstrate T
if sSimParams.b_plotTransDemos
    for d = 1:dim
        PlotPolyCdfDemonstration1(xMin(d), xMax(d), pCdf(d,:), xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), muTilde(d), sigmaTilde(d,d));
        PlotPolyCdfDemonstration2(xMin(d), xMax(d), pCdf(d,:), invpCdf(d,:), muTilde(d), sigmaTilde(d,d));
    end
    if dim == 2
        PlotPolyCdfDemonstration3_2D(xTrain,xTildeTrain)
    end
end
%% G tilde parameters
omegaTilde = sDataset.recommendedOmega;
MTilde = 50;

%% Calculate analytic eigenfunctions of W_tilde
[PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
sqrtnLambdaPhiTilde = sqrt(n*lambdaAnalyticTilde)'.*PhiTilde;

if sSimParams.b_plotAllEvecs
    % Build W tilde and numeric eigenvectors for comparison
    distTilde = pdist2(xTildeTrain, xTildeTrain);
    WTilde = exp(-distTilde.^2/(2*omegaTilde^2));
    [VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
    lambdaNumericTilde = diag(LambdaNumericTilde);
    sqrtLambdaVTilde = abs(sqrt(lambdaNumericTilde))'.*VTilde;
    VTilde = FlipSign(PhiTilde, VTilde);

    figTitle = 'Eigenfunctions \& eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    figName = 'PhiTilde_VTilde';
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1), vInd(end), ...
        sqrtnLambdaPhiTilde, lambdaAnalyticTilde, [], [], figTitle, figName, ...
        '\sqrt{n \tilde{\lambda}^{\phi}}\tilde{\phi}', sqrtLambdaVTilde, '\sqrt{\tilde{\lambda}^{v}}\tilde{v}')
    
    PlotSpectrum(sSimParams, sDataset, [], n*lambdaAnalyticTilde, lambdaNumericTilde, [], 'n \tilde{\lambda}^{\phi}_m', '\tilde{\lambda}^{v}_m');

    vPrTilde = SimpleEstPorbablityArea(xTildeTrain, sigmaTilde, muTilde);
    PlotInnerProductMatrix(sSimParams, dim, vPrTilde, 'IP_Matrix', [], PhiTilde, 'Analytic');
    PlotInnerProductMatrix(sSimParams, dim, [], 'IP_Matrix', [], PhiTilde, 'Numeric');
end
%% Functional maps
C = pinv(PhiTilde)*V;

figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

%% V in terms of Phi_tilde
if sSimParams.b_plotVRec
    VRec = PhiTilde*C;

    figTitle = '${\bf V}$ in terms of ${\bf \tilde{\Phi}} \quad ({\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C})$';
    figName = 'VRec';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1), vInd(end), ...
        VRec, lambda, [], [], figTitle, figName, 'v^{{\bf rec}}')

    b_plotErrVsNodeInd = true;
    PlotEigenDiffs(sSimParams, sDataset, [], vInd(1), vInd(end), VRec, V, 'VRec_vs_V', 'v^{{\bf rec}}', 'v', b_plotErrVsNodeInd)
end
%% Interpolate
xInt = sDataset.sData.xt;

xTildeInt = T(pCdf, true, muTilde, sigmaTilde, xInt);
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);
% PhiTildeInt = sqrt(n*lambdaAnalyticTilde)'.*PhiTildeInt;
if sSimParams.b_plotAllEvecs
    figName = 'PhiTildeInt';
    figTitle = GetInterpFigTitle(dim, omegaTilde, sigmaTilde, muTilde);
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeInt, [], vInd(1), vInd(end), ...
        PhiTildeInt, lambdaAnalyticTilde, [], [], figTitle, figName, '\tilde{\phi}^{{\bf int}}')
end
% xIntInvT = invT(invpCdf, muTilde, sigmaTilde, xTildeInt);
interpRatio = N/n;
VInt = sqrt(interpRatio)*PhiTildeInt*C;

if sSimParams.b_plotOrigVsInterpEvecs || sSimParams.b_plotAllEvecs
    figTitle = 'Interpolated eigenvectors of ${\bf W}$ - Ours';
    figName = 'VInt';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
        VInt, lambda, [], [], figTitle, figName, 'v^{{\bf int}}')
end
%% Interpolate with Nystrom
if strcmp(interpMethod, 'AddPoints')
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omega^2));

    VNys = B.'*V*diag(1./lambda);
%     VNys = (1/sqrt(interpRatio))*VNys;
    if sSimParams.b_plotOrigVsInterpEvecs || sSimParams.b_plotAllEvecs

        figTitle = 'Interpolated eigenvectors of ${\bf W}$ - Nystrom';
        figName = 'VNys';
        PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
            VNys, lambda, [], [], figTitle, figName, 'v^{{\bf Nys}}')
    end
else
    warning('Cannot run Nystrom when interpMethod is not ''AddPoints''')
end

%% Reference
WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega, k);
[VRef, LambdaRef] = eigs(WRef,M);
lambdaRef = diag(LambdaRef);
VRef = FlipSign(VInt, VRef);
% VRef = (sqrt(lambdaRef))'.*VRef;

if sSimParams.b_plotOrigVsInterpEvecs || sSimParams.b_plotAllEvecs
    figTitle = 'Reference eigenvectors of ${\bf W}^{{\bf ref}}$';
    figName = 'VRef';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
        VRef, lambdaRef, [], [], figTitle, figName, 'v^{{\bf ref}}')
end
%% Compare
lambdaRenormedToCompare = interpRatio*lambda;
PlotSpectrum(sSimParams, sDataset, [], lambdaRenormedToCompare, lambdaRef, [], '\frac{N}{n} \lambda_m', '\lambda^{{\bf ref}}');

VRenormedToCompare = (sqrt(lambda))'.*V;
VRefRenormedToCompare = (sqrt(lambdaRef))'.*VRef;
figTitle = '$\sqrt{\Lambda} V$ vs. $\sqrt{\Lambda^{{\bf ref}}} V^{{\bf ref}}_{(1:n,:)}$';
figName = 'V_vs_VRef';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1), vInd(end), ...
    VRenormedToCompare, [], [], [], figTitle, figName, '\sqrt{\lambda} v', VRefRenormedToCompare(1:n,:), '\sqrt{\lambda}^{{\bf ref}} v^{{\bf ref}}');

VIntRenormedToCompare = (sqrt(lambda))'.*VInt;
figTitle = 'Ours vs. Reference';
figName = 'VInt_vs_VRef';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VIntRenormedToCompare, lambdaRef, [], [], figTitle, figName, 'v^{{\bf int}}', VRefRenormedToCompare, 'v^{{\bf ref}}');


VNysRenormedToCompare = (sqrt(lambda))'.*VNys;
figTitle = 'Nystrom vs. Reference';
figName = 'VNys_vs_VRef';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VNysRenormedToCompare, lambdaRef, [], [], figTitle, figName, 'v^{{\bf nys}}', VRefRenormedToCompare, 'v^{{\bf ref}}');


figTitle = 'Ours vs. Nystrom';
figName = 'VInt_vs_VNys';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VInt, lambdaRef, [], [], figTitle, figName, 'v^{{\bf int}}', VNys, 'v^{{\bf nys}}');

figTitle = 'Ours vs. Nystrom';
figName = 'VIntRenormedToCompare_vs_VNysRenormedToCompare';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VIntRenormedToCompare, lambdaRef, [], [], figTitle, figName, 'v^{{\bf int}}', VNysRenormedToCompare, 'v^{{\bf nys}}');


b_plotErrVsNodeInd = true;
PlotEigenDiffs(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), VIntRenormedToCompare, VRefRenormedToCompare, ...
    'VInt_vs_VRef', 'v^{{\bf int}}', 'v^{{\bf ref}}', b_plotErrVsNodeInd)
PlotEigenDiffs(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), VNysRenormedToCompare, VRefRenormedToCompare, ...
    'VNys_vs_VRef', 'v^{{\bf nys}}', 'v^{{\bf ref}}', b_plotErrVsNodeInd)

%% RMSE
R = 10;
mVIntToCompare = zeros(R, N, M);
mVNysToCompare = zeros(R, N, M);
mVRefToCompare = zeros(R, N, M);
for r = 1:R
    % Generate dataset
    sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod);
    xTrain = sDataset.sData.x;
    xInt = sDataset.sData.xt;
    interpRatio = N/n;
    
    % Original graph
    W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega, k);
    [V, Lambda] = eigs(W,M);
    
    % Transform with pCdf, calc eigen functions and C
    [~, ~, pCdf, invpCdf] = PolyfitEstCdf(xTrain, bins, pCdfDegree, invpCdfDegree, b_kde);
    xTildeTrain = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain);
    [PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
    C = pinv(PhiTilde)*V;
    
    % Interpolate with our method
    xTildeInt = T(pCdf, true, muTilde, sigmaTilde, xInt);
    [PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);
    VInt = sqrt(interpRatio)*PhiTildeInt*C;
%     VIntToCompare = abs((sqrt(lambda))'.*VInt);
    VIntToCompare = abs(VInt);
    
    % Interpolate with Nystrom
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omega^2));
    VNys = B.'*V*diag(1./lambda);
%     VNysToCompare = abs((sqrt(lambda))'.*VNys);
    VNysToCompare = abs(VNys);

    % Calculate reference
    WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega, k);
    [VRef, LambdaRef] = eigs(WRef,M);
    VRef = FlipSign(VInt, VRef);
%     VRefToCompare = abs(abs((sqrt(lambdaRef)))'.*VRef);
    VRefToCompare = sqrt(interpRatio)*abs(VRef);
    if r == 1
        figTitle = 'Ours vs. Reference';
        figName = 'VInt_vs_VRef';
        PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], 22, 24, ...
            VIntToCompare, lambdaRef, [], [], figTitle, figName, 'v^{{\bf int}}', VRefToCompare, 'v^{{\bf ref}}');

        figTitle = 'Nystrom vs. Reference';
        figName = 'VNys_vs_VRef';
        PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], 22, 24, ...
            VNysToCompare, lambdaRef, [], [], figTitle, figName, 'v^{{\bf nys}}', VRefToCompare, 'v^{{\bf ref}}');
    end
    
    % Save
    mVIntToCompare(r,:,:) = VIntToCompare;
    mVNysToCompare(r,:,:) = VNysToCompare;
    mVRefToCompare(r,:,:) = VRefToCompare;

end
vRmseInt = CalcRMSE(mVIntToCompare, mVRefToCompare, 'Analytic');
vRmseNys = CalcRMSE(mVNysToCompare, mVRefToCompare, 'Nystrom');

windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')
figure;
% subplot(2,1,1)
plot((0:M-1)', [vRmseInt.' vRmseNys.'], 'LineWidth', 2)
ylim([0 0.5]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Ours', 'Nystrom', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);
% subplot(2,1,2)
% plot((0:M-1)', [10*log10(vRmseInt).' 10*log10(vRmseNys).'], 'LineWidth', 2)
% % set(gca, 'YScale', 'log')
% xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
% legend('Ours', 'Nystrom', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
% set(gca,'FontSize', 14);
set(0,'DefaultFigureWindowStyle',windowStyle)
