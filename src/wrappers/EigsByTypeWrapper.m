function [V, VRef, adjLambda, matLambda, adjLambdaRef, matLambdaRef] = EigsByTypeWrapper(sPlotParams, sPreset, sDataset, W, D, Ln, WRef, DRef, LnRef)
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
if isfield(sDataset.sData, 'S')
    S = sDataset.sData.S;
    SInt = sDataset.sData.St;
else
    S = sDataset.sData.x;
    SInt = sDataset.sData.xt;
end
if b_debugUseAnalytic
    fprintf('Generating analytic expression for %s %s\n', verticesPDF, matrixForEigs);
    if ismember(verticesPDF, {'Grid', 'Uniform', 'SwissRoll'})
        assert(strcmp(matrixForEigs, 'Laplacian') || strcmp(matrixForEigs, 'RandomWalk'))
        if ismember(verticesPDF, {'Grid', 'Uniform'})
            len = sDatasetParams.xMax - sDatasetParams.xMin;
        else % Swiss Roll
            len = [SwissRollArclength(sDatasetParams.maxTheta), sDatasetParams.height];
        end
        [V, matLambda] = CalcAnalyticLapEigsGrid(S, M, len);
        [VRef, matLambdaRef] = CalcAnalyticLapEigsGrid(SInt, M, len);
        [~, adjLambda, ~] = EigsByType(W, D, Ln, M, matrixForEigs);
        [~, adjLambdaRef, ~] = EigsByType(WRef, DRef, LnRef, M, matrixForEigs);
    elseif strcmp(verticesPDF, 'Gaussian')
        assert(strcmp(matrixForEigs, 'Adjacency') && dim == 1)
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
if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotOrigEvecs && dim <= 3
    if dim == 1
        plotInd = [0,min(4,M-1)];
    else
        plotInd = [0,min(8,M-1)];
    end
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
    if isfield(sDataset.sData, 'S')
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, S, [], plotInd(1), plotInd(end), ...
            V, matLambda, '\lambda^{v}', [], 'Parametrization', figName, 'v');
    end
    figTitle = 'Numeric eigenvalues of ${\bf V}$';
    PlotSpectrum([], [], matLambda, [], [], '\tilde{\lambda}^{v}_m', [], [], figTitle);
end
end