function [] = PlotGaussianKernelEigenfunsExample()
sigma = 1;
omega = sqrt(sigma/10); 0.3;
mu = 5;

n = 2000;
firstEigInd = 0;
lastEigInd = 5;
M = lastEigInd+1;

x = mu + linspace(-3*sigma, 3*sigma, n).';
sPreset.dim = 1;
sPreset.verticesPDF = 'Gaussian';
sPreset.matrixForEigs = 'Adjacency';
b_saveFigures = true;
sPlotParams = GetPlotParams(sPreset, b_saveFigures);

sDistParams.dim = 1;
sDistParams.cov{1} = sigma;
sDistParams.sigma{1} = diag(sigma);
sDistParams.sigma_eigv{1} = sigma;
[sDistParams.u{1}, sDistParams.sigma_eigv{1}] = eig(sDistParams.cov{1});
sDistParams.mu{1} = mu;
sDistParams.mu_1D{1} = sDistParams.mu{1}*sDistParams.u{1};
sDistParams.estNumComponents = 1;
sDistParams.componentProportion = 1;
sDistParams.estDataDist = 'Gaussian';
sDistParams.GMModel.ComponentProportion = 1;

sKernelParams = CalcKernelParams(sDistParams, omega);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(M, sKernelParams);
[ Phi, lambda ] = CalcAnalyticEigenfunctions(M, sKernelParams, x);

figTitle = 'Eigenfunctions of the Gaussian kernel on $\mathbb{R}$';
plotInd = firstEigInd:lastEigInd;
[cData{1:numel(plotInd)}] = deal(x);
cSigStr = RepLegend('\\phi', plotInd);
[cNumCircles{1:numel(plotInd)}] = deal((1:n).');
[cMarkers{1:numel(plotInd)}] = deal('.');
PlotGraphSignals(sPlotParams, [], ...
    ['PhiTilde_' ,num2str(plotInd(1)), '_to_', num2str(plotInd(end)) ], cData, ...
    [mat2cell(Phi(:,plotInd+1),n,ones(1,numel(plotInd))), mat2cell(Phi(:,plotInd+1),n,ones(1,numel(plotInd)))], ...
    cSigStr, cNumCircles, cMarkers, [], [], [min(min(Phi(:,plotInd+1))), max(max(Phi(:,plotInd+1)))]);

% figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
PlotSpectrum([], [], lambda, [], [], '\lambda_m', [], [], []);
end