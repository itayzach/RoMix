function [] = PlotGaussianKernelEigenfunsExample()
sigma = 1;
omega = sqrt(sigma/10); 0.3;
mu = 5;

n = 2000;
firsEigInd = 0;
lastEigInd = 5;
M = lastEigInd+1;

b_normalizePhi = false;
x = mu + linspace(-3*sigma, 3*sigma, n).';
sPreset.dim = 1;
sPreset.verticesPDF = 'Gaussian';
sPreset.matrixForEigs = 'Adjacency';
sPlotParams = GetPlotParams(sPreset);

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

sKernelParams = CalcKernelParams(sDistParams, omega);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(M, sKernelParams);
[ Phi, lambda ] = CalcAnalyticEigenfunctions(M, sKernelParams, x, b_normalizePhi);

figTitle = 'Eigenfunctions of the Gaussian kernel on $\mathbb{R}$';
figName = 'PhiTilde';
PlotEigenfuncvecScatter(sPlotParams, 'Gaussian', x, [], firsEigInd, lastEigInd, ...
    Phi, [], [], [], [], figName, '\phi' );
% figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
PlotSpectrum([], [], lambda, [], [], '\lambda_m', [], [], []);
end