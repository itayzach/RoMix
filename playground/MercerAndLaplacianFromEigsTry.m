clc; clear; close all;

MTilde = 80;
omega = 0.07;
sPreset = Get1DGridPreset();
sDataset = GenerateDataset(sPreset.verticesPDF, sPreset.dim, sPreset.nGenDataCompnts, ...
    sPreset.n, sPreset.N, 'AddPoints', sPreset.sDatasetParams);

%% Mercer
sDistParams = EstimateDistributionParameters(sDataset.sData.x, ...
    sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);
sKernelParams = GetKernelParams(sDistParams, omega);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(MTilde, sKernelParams);
[ PhiTilde, lambdaAnalyticTilde ] = ...
    CalcAnalyticEigenfunctions(MTilde, sKernelParams, sDataset.sData.x, sPreset.b_normalizePhi);
WAnalyticMercer = (PhiTilde*diag(lambdaAnalyticTilde)*PhiTilde');

%% Numeric
[W, dist, D] = SimpleCalcAdjacency(sDataset.sData.x, sPreset.adjacencyType, ...
    sPreset.distType, omega, sPreset.k, sPreset.nnValue);

%% Plot
% sPlotParams = GetPlotParams();
% PlotEigenfuncvecScatter(sPlotParams, 'Gaussian', sDataset.sData.x, [], MTilde-1, MTilde-1, ...
%     PhiTilde, [], [], [], 'Eigenfunctions', 'efuncs', '\tilde{\phi}' );

% PhiTildeFft = fftshift(fft(PhiTilde,[],1));
% figure;
% plot([abs(PhiTildeFft(:,1)), abs(PhiTildeFft(:,end))])

figure;
plot([PhiTilde(:,1), PhiTilde(:,end)])

figure;
tiledlayout('flow');
ax(1) = nexttile();
    imagesc(WAnalyticMercer);
    colorbar()
    caxis([0 max(1,max(WAnalyticMercer(:)))])
    title('${\bf W}_{i,j} = \sum_{m=1}^M \lambda_m \phi_m(x_i)\phi_m(x_j), \quad {\bf W} = \Phi \Lambda \Phi^T$', ...
        'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
ax(2) = nexttile();
    imagesc(W);
    colorbar()
    caxis([0 1])
    title('${\bf W}_{i,j} = \exp\big(\|x_i-x_j\|^2/\omega^2\big)$', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
linkaxes(ax,'xy')
%% Verify
threshold = 1e-3;
b_AssertOnFail = true;
b_printMaxErr = true;
isalmostequal(W, WAnalyticMercer, threshold, [], b_AssertOnFail, b_printMaxErr)




%% Laplacian try
L = D - W;
% d = diag(D);
% Wn = diag(d.^-0.5)*W*diag(d.^-0.5);
% I = eye(length(W));
% Ln = I - Wn;
% [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, DRef, M, matrixForEigs);
LAnalytic = (PhiTilde*diag(exp(-sKernelParams.t*flipud(lambdaAnalyticTilde)))*PhiTilde');

figure;
tiledlayout('flow');
ax(1) = nexttile();
    imagesc(LAnalytic);
    colorbar()
    title('LAnalytic', ...
        'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
ax(2) = nexttile();
    imagesc(L);
    colorbar()
    title('L', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
linkaxes(ax,'xy')

