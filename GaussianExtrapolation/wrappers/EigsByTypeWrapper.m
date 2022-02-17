function [V, adjLambda, matLambda, VRef, adjLambdaRef, matLambdaRef] = EigsByTypeWrapper(sPlotParams, sPreset, sDataset, W, D, Ln, WRef, DRef, LnRef)
dim                = sPreset.dim;
n                  = sPreset.n;
N                  = sPreset.N;
k                  = sPreset.k;
verticesPDF        = sPreset.verticesPDF;
adjacencyType      = sPreset.adjacencyType;
matrixForEigs      = sPreset.matrixForEigs;
M                  = sPreset.M;
omega              = sPreset.omega; % for W
sDatasetParams     = sPreset.sDatasetParams;
b_debugUseAnalytic = sPreset.b_debugUseAnalytic;
b_takeEigsFromWRef = sPreset.b_takeEigsFromWRef;
interpRatio        = N/n;

xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
%%
if b_debugUseAnalytic
    assert(strcmp(matrixForEigs, 'Adjacency') && dim == 1)
    if strcmp(verticesPDF, 'Grid')
        % analytic eigenfunctions from Spectral Graph Theory, Spielman
        V = (1/sqrt(n))*cos(pi*(0:M-1).*xTrain(:,1) - pi*(0:M-1));
        VRef = (1/sqrt(N))*cos(pi*(0:M-1).*xInt(:,1) - pi*(0:M-1));
        adjLambda = 2*(1-cos(pi*(0:M-1)/n));
        adjLambdaRef = 2*(1-cos(pi*(0:M-1)/N));
        matLambda = adjLambda;
        matLambdaRef = adjLambdaRef;
    elseif strcmp(verticesPDF, 'Gaussian')
        warning('review the following {1}...')
        [V, adjLambda] = SimpleCalcAnalyticEigenfunctions(xTrain, omega, sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        [VRef, adjLambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        adjLambdaRef = N*adjLambdaRef;
        matLambda = adjLambdaRef;
    else
        error('invalid pdf option for analytic expressions')
    end
else
    [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, DRef, LnRef, M, matrixForEigs);
    if b_takeEigsFromWRef
        V = sqrt(interpRatio)*VRef(1:n,:);
        adjLambda = adjLambdaRef/interpRatio;
        matLambda = matLambdaRef/interpRatio;
    else
        [V, adjLambda, matLambda] = EigsByType(W, D, Ln, M, matrixForEigs);
    end
end
assert(isreal(V), 'V should be real...')
assert(isreal(matLambda), 'matLambda should be real...')

% ----------------------------------------------------------------------------------------------
% Plot numeric V
% ----------------------------------------------------------------------------------------------
if dim == 1
    plotInd = [0,min(4,M-1)];
else
    plotInd = [0,min(8,M-1)];
end
if sPlotParams.b_globalPlotEnable && (sPlotParams.b_plotOrigEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
    figTitle = ['Eigenvectors of ', matrixForEigs];
    if strcmp(adjacencyType, 'GaussianKernel')
        figTitle2 = [' (Gaussian kernel, $\omega = ', num2str(omega), '$'];
    elseif strcmp(adjacencyType, 'NearestNeighbor')
        figTitle2 = [' (k-NN, $k = ', num2str(k), '$'];
    end
    figTitle3 = [', $n = ', num2str(n), '$)'];
    figName = 'V';
    PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
        V, matLambda, '\lambda^{v}', [], [figTitle, figTitle2, figTitle3], figName, 'v');
    figTitle = 'Numeric eigenvalues of ${\bf V}$';
    PlotSpectrum([], [], matLambda, [], [], '\tilde{\lambda}^{v}_m', [], [], figTitle);
end
end