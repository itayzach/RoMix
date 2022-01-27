%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')

%%
sPlotParams = [];
sPreset = Get1DGridPreset();
%%
dim                = sPreset.dim;
n                  = sPreset.n;
N                  = sPreset.N;
k                  = sPreset.k;
nnValue            = sPreset.nnValue;
verticesPDF        = sPreset.verticesPDF;
adjacencyType      = sPreset.adjacencyType;
matrixForEigs      = sPreset.matrixForEigs;
nGenDataCompnts    = sPreset.nGenDataCompnts;
dataGenTechnique   = sPreset.dataGenTechnique;
M                  = sPreset.M;
MTilde             = sPreset.MTilde;
omega              = sPreset.omega;
omegaNys           = sPreset.omega;
omegaRep           = sPreset.omega;
omegaTilde         = sPreset.omegaTilde;
gamma1             = sPreset.gamma1;
gamma2             = sPreset.gamma2;
gamma1Rep          = sPreset.gamma1Rep;
gmmRegVal          = sPreset.gmmRegVal;
gmmMaxIter         = sPreset.gmmMaxIter;
gmmNumComponents   = sPreset.gmmNumComponents;
sDatasetParams     = sPreset.sDatasetParams;
sDistanceParams    = sPreset.sDistanceParams;
b_normalizePhi     = sPreset.b_normalizePhi;
b_maskDataFitTerm  = false;
plotInd = [0,min(5,M-1)];
%%
sDataset = GenerateDataset(verticesPDF, dim, nGenDataCompnts, n, N, dataGenTechnique, sDatasetParams);
%%
x = sort(sDataset.sData.x);
xInt = sort(sDataset.sData.xt);
%% Signal
fs = N/(sDatasetParams.xMax-sDatasetParams.xMin);
f0 = 1;
signal = sin(2*pi*f0*x);
signalRef = sin(2*pi*f0*xInt);

%% Classic DFT
MDft = 20;
tic;
VDftComplex = (1/sqrt(N))*dftmtx(N);
t = toc;
VDft = VDftComplex(:,1:MDft);
fprintf('dftmtx took %.2f sec\n',t);
% VDft(:,1:2:2*MDft) = real(VDftComplex(:,1:MDft));
% VDft(:,2:2:2*MDft) = imag(VDftComplex(:,1:MDft));
dftSignal = VDft.'*signalRef; % <v_k,signal>_{\ell_2} --> p(x)=uniform
% freqs = (0:2*MDft-1)/N * fs;
freqs = (0:MDft-1)/N * fs;

figure;
subplot(3,1,1)
plot(xInt,signalRef,'o'); hold on; plot(x,signal,'.');
title('$s(x)$', 'interpreter', 'latex');  set(gca,'FontSize', 14);
xlabel('x [sec]', 'interpreter', 'latex')
set(gca,'FontSize', 14);
subplot(3,1,2);
% plot(xInt, VDft(:,1:6), '.')
plot([real(VDft(:,1:3)), imag(VDft(:,1:3))], '.')
title('DFT matrix, $v_k(x)$', 'interpreter', 'latex'); set(gca,'FontSize', 14);
subplot(3,1,3)
stem(freqs, abs(dftSignal), 'filled')
xlabel('f [Hz]', 'interpreter', 'latex')
title('$|$DFT$\{s\}|=|{\bf V}_{{\bf DFT}}^T s|$', 'interpreter', 'latex'); set(gca,'FontSize', 14);
% x0 = 100; y0 = 100; width = 600; height = 600;
% set(gcf,'Position', [x0 y0 width height])

%% Build graph
% [WRef, distRef, DRef, DRefsqinv] = SimpleCalcAdjacency(xInt, adjacencyType, sDistanceParams, omega, k, nnValue);
% [W, dist, D, Dsqinv] = SimpleCalcAdjacency(x, adjacencyType, sDistanceParams, omega, k, nnValue);
%% Make sure x and xInt are sorted
% WRef = toeplitz([ 0 1 zeros(1, N-2)]); % path graph
% WRef = toeplitz([ 1 1 zeros(1, N-2)]); % path graph with self loops
% WRef = toeplitz([ 0 1 zeros(1, N-3) 1]); % ring graph
WRef = toeplitz([ 1 1 zeros(1, N-3) 1]); % ring graph with self loops

W = toeplitz([ 1 1 zeros(1, n-3) 1]); % ring graph with self loops
%% Laplacian
dRef = sum(WRef);
DRefsqinv = diag(1./sqrt(dRef));
LnRef = eye(N) - DRefsqinv*WRef*DRefsqinv;
PlotWeightsMatrix([],WRef,[],dRef,LnRef,xInt,adjacencyType,omega,k)

d = sum(W);
Dsqinv = diag(1./sqrt(d));
Ln = eye(n) - Dsqinv*W*Dsqinv;
PlotWeightsMatrix([],W,[],d,Ln,x,adjacencyType,omega,k)

%% GFT with eigs(Ln)
MGft = 20;
tic;
[URefGsp, LambdaGsp] = eig(LnRef);
t = toc;
fprintf('eig(LnRef) took %.2f sec\n',t);
URefGsp = URefGsp(:,1:MGft);
LambdaGsp = LambdaGsp(1:MGft,1:MGft);
gftSignalRef = URefGsp.'*signalRef;
lambdaGsp = diag(LambdaGsp);

PlotGraphSignals(sPlotParams, '$s(x)$', 'signal', {x, xInt}, {signal, signalRef}, ...
    {'$s(x)$', '$s^{{\bf ref}}(x)$'});
PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
    URefGsp, [], [], [], ...
    ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf ref}}');
figure;
stem(abs(gftSignalRef), 'filled')
xlabel('k', 'interpreter', 'latex')
title('GFT$\{s\} = |{\bf U}_{{\bf GFT}}^T s|$', 'interpreter', 'latex'); set(gca,'FontSize', 14);
% x0 = 700; y0 = 100; width = 600; height = 200;
% set(gcf,'Position', [x0 y0 width height])

%% GFT with our approach
tic;
[UGsp, ~] = eig(Ln);
UGsp = UGsp(:,1:MGft);
sDistParams = EstimateDistributionParameters(x, gmmNumComponents, gmmRegVal, gmmMaxIter);
sKernelParams = CalcKernelParams(sDistParams, omega);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(MTilde, sKernelParams);
[ Phi, lambda ] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, x, b_normalizePhi);
invLambda = diag(1./lambda);
C = EigsRLS(Phi, gamma1, 0, invLambda, [], UGsp, b_maskDataFitTerm);
C = C/sqrt(N/n);
[ PhiInt, ~ ] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xInt, b_normalizePhi);
t = toc;
fprintf('our approach took %.2f sec\n',t);
UOurs = PhiInt*C;
UOurs = FlipSign(URefGsp,UOurs);
oursGftSignalInt = UOurs.'*signalRef;

PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
    UOurs, [], [], [], ...
    ['Ours (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf int}}');
figure;
stem(abs(oursGftSignalInt), 'filled')
xlabel('k', 'interpreter', 'latex')
title('OurGFT$\{s\} = |({\bf \Phi}{\bf C})^T s|$', 'interpreter', 'latex'); set(gca,'FontSize', 14);
% x0 = 1300; y0 = 100; width = 600; height = 600;
% set(gcf,'Position', [x0 y0 width height])

% vRmseNys = CalcRMSE(mVNysToCompare, oursGftSignal, 'Nystrom');
% vRmseRep = CalcRMSE(mVRepToCompare, mVRefToCompare, 'Representer');
% PlotRMSE(sPlotParams, [vRmseInt, vRmseNys, vRmseRep], ...
%     {'RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
%     'RMSE$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'});
