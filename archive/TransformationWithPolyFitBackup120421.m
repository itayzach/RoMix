%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')
sSimParams.outputFolder               = 'figs';
sSimParams.b_plotAllEvecs             = true;
sSimParams.b_GSPBoxPlots              = false;
sSimParams.b_plotLaplacianEvecs       = false;
sSimParams.b_plotWeights              = false;
sSimParams.b_plotVRec                 = false;
sSimParams.b_plotTransDemos           = true;
sSimParams.b_plotOrigVsInterpEvecs    = true;

%% Random data
dim = 1;
nComponents = 1;
n = 2000;
N = 4000;
verticesPDF = 'Uniform'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
origGraphAdjacency = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
interpMethod = 'AddPoints'; % 'NewPoints' / 'AddPoints'
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod);
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
k = 3;
W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega, k);

%% Calculate (numeric) eigenvectors of G
M = 30;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
% vInd = 0:4;
vInd = M-5:M-1;

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
[xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = PolyfitEstCdf(xTrain, bins, pCdfDegree, invpCdfDegree);

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
%% Build G tilde and numeric eigenvectors for comparison
omegaTilde = sDataset.recommendedOmega;
MTilde = 30;
distTilde = pdist2(xTildeTrain, xTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);
%% Calculate analytic eigenfunctions of W_tilde
[PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);

if sSimParams.b_plotAllEvecs
    figTitle = 'Eigenfunctions \& eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
    figName = 'PhiTilde_VTilde';
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1), vInd(end), ...
        PhiTilde, lambdaAnalyticTilde, [], [], figTitle, figName, '\tilde{\phi}', VTilde, '\tilde{v}')
    
    vPrTilde = SimpleEstPorbablityArea(xTildeTrain, sigmaTilde, muTilde);
    PlotInnerProductMatrix(sSimParams, dim, vPrTilde, 'IP_Matrix', [], PhiTilde, 'Analytic');
    PlotInnerProductMatrix(sSimParams, dim, [], 'IP_Matrix', [], PhiTilde, 'Numeric');
    
    % STAM
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeTrain, [], vInd(1), vInd(end), ...
        V, lambdaAnalyticTilde, [], [], figTitle, figName, '\tilde{\phi}', VTilde, '\tilde{v}')
    
    
end
%% Functional maps
% C = diag(sqrt(lambdaAnalyticTilde))*diag(sqrt(lambda))*pinv(PhiTilde)*V;
C = pinv(PhiTilde)*V;
% C = PhiTilde.'*V;
% error('take sqrt(lambda) into account?')

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

if sSimParams.b_plotAllEvecs
    figName = 'PhiTildeInt';
    figTitle = GetInterpFigTitle(dim, omegaTilde, sigmaTilde, muTilde);
    PlotEigenfuncvecScatter(sSimParams, 'Gaussian', xTildeInt, [], vInd(1), vInd(end), ...
        PhiTildeInt, lambdaAnalyticTilde, [], [], figTitle, figName, '\tilde{\phi}^{{\bf int}}')
end
% xIntInvT = invT(invpCdf, muTilde, sigmaTilde, xTildeInt);


interpRatio = N/n;
VIntRenormed = sqrt(interpRatio)*PhiTildeInt*C;

if sSimParams.b_plotOrigVsInterpEvecs || sSimParams.b_plotAllEvecs
    figTitle = 'Interpolated eigenvectors of ${\bf W}$ - Ours';
    figName = 'VInt';
    PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
        VIntRenormed, lambda, [], [], figTitle, figName, 'v^{{\bf int}}')
end
%% Interpolate with Nystrom
if strcmp(interpMethod, 'AddPoints')
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omega^2));

    VNys = B.'*V*diag(1./lambda);
    VNys = (1/sqrt(interpRatio))*VNys;
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
VRefRenormedToCompare = (sqrt(lambdaRef))'.*VRef(1:n,:);
figTitle = '$\sqrt{\Lambda} V$ vs. $\sqrt{\Lambda^{{\bf ref}}} V^{{\bf ref}}$ on n points';
figName = 'V_vs_VRef';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xTrain, [], vInd(1), vInd(end), ...
    VRenormedToCompare, [], [], [], figTitle, figName, '\sqrt{\lambda} v', VRefRenormedToCompare, '\sqrt{\lambda}^{{\bf ref}} v^{{\bf ref}}');

figTitle = 'Ours vs. Nystrom';
figName = 'VInt_vs_VNys';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VInt, lambdaRef, [], [], figTitle, figName, 'v^{{\bf int}}', VNys, 'v^{{\bf nys}}');

figTitle = 'Ours vs. Reference';
figName = 'VInt_vs_VRef';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VInt, lambdaRef, [], [], figTitle, figName, 'v^{{\bf int}}', VRef, 'v^{{\bf ref}}');

figTitle = 'Nystrom vs. Reference';
figName = 'VNys_vs_VRef';
PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), ...
    VNys, lambdaRef, [], [], figTitle, figName, 'v^{{\bf nys}}', VRef, 'v^{{\bf ref}}');

b_plotErrVsNodeInd = true;
PlotEigenDiffs(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), VInt, VRef, ...
    'VInt_vs_VRef', 'v^{{\bf int}}', 'v^{{\bf ref}}', b_plotErrVsNodeInd)
PlotEigenDiffs(sSimParams, verticesPDF, xInt, [], vInd(1), vInd(end), VNys, VRef, ...
    'VNys_vs_VRef', 'v^{{\bf nys}}', 'v^{{\bf ref}}', b_plotErrVsNodeInd)

%% RMSE
R = 10;
mVInt = zeros(R, N, M);
mVNys = zeros(R, N, M);
mVRef = zeros(R, N, M);
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
    [~, ~, pCdf, invpCdf] = PolyfitEstCdf(xTrain, bins, pCdfDegree, invpCdfDegree);
    xTildeTrain = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain);
    
    % Interpolate with our method
    [PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
    C = pinv(PhiTilde)*V;
    VInt = PhiTildeInt*C;
    
    % Interpolate with Nystrom
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omega^2));
    VNys = B.'*V*diag(1./lambda);
    VNys = (1/sqrt(interpRatio))*VNys;
    
    % Calculate reference
    WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega, k);
    [VRef, LambdaRef] = eigs(WRef,M);
    VRef = FlipSign(VInt, VRef);
    
    % Save
    mVInt(r,:,:) = VInt;
    mVNys(r,:,:) = VNys;
    mVRef(r,:,:) = VRef;

end
vRmseInt = CalcRMSE(mVInt, mVRef, 'Analytic');
vRmseNys = CalcRMSE(mVNys, mVRef, 'Nystrom');

windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')
figure;
plot((0:M-1)', [vRmseInt.' vRmseNys.'], 'LineWidth', 2)
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Ours', 'Nystrom', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);
set(0,'DefaultFigureWindowStyle',windowStyle)
