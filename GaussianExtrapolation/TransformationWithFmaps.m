%% Restart
clc; clear; close all; rng('default')

%% Parameters
%==========================================================================
% Graph data
%==========================================================================
% 'Gaussian_1D' / 'Uniform_1D' / 'Gaussian_2D' / 'Uniform_2D' / 
% 'Two_moons' /
% 'bunny' / 'twodgrid' / 'sensor' / 'minnesota' /
% 'david_sensor' / 'swiss_roll' / 'random_ring'
graphName =  'Uniform_2D'; 

%==========================================================================
% Graph & graph-signal parameters
%==========================================================================
omega       = 0.3;
estDataDist = 'Gaussian';
nComponents = 1;
b_normlizedLaplacian = true;
graphSignalModel = 'V_c'; % 'alpha_K' / 'V_c' / 'Phi_c'
G_tildeBasis = 'Numeric'; % 'Analytic' / 'Numeric'

nEigs = 30;
samplingRatio = 0;
nysRatio = 0.8;

%==========================================================================
% Plot parameters
%==========================================================================
if samplingRatio < 0.2
    sSimParams.b_plotSamplingPointsMarkers = true;
else
    sSimParams.b_plotSamplingPointsMarkers = false;
end
sSimParams.outputFolder               = 'figs';
sSimParams.PlotEigenFuncsM            = min(nEigs, 12);
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = true;
%% Generate graph
[G, dist, sDataset, sKernelParams] = GenerateGraph(graphName, nComponents, ...
    estDataDist, omega, b_normlizedLaplacian, nEigs);
v = sDataset.sData.x; % For convenience
if sSimParams.b_showGraphMatricesFigures
    PlotGraphMatrices(G, b_normlizedLaplacian);
end
if sDataset.dim == 1
    sSimParams.b_GSPBoxPlots = false;
end

%==========================================================================
% Eigenvectors
%==========================================================================
[V, vLambdaNumeric] = CalcNumericEigenvectors(nEigs, sKernelParams, v);
if sSimParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, ...
        sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [], G, ...
        '${\bf W}_G$ Eigenvectors')
end

% %==========================================================================
% % Nystrom Eigenvectors
% %==========================================================================
% [Vnys, vLambdaNys] = CalcNystromEigenvectors(nEigs, sKernelParams, v, nysRatio);
% if sSimParams.b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v, nysRatio, 0, ...
%         sSimParams.PlotEigenFuncsM-1, Vnys, vLambdaNys, 'Nystrom', [], G, ...
%         '${\bf W}_G$ Nystrom Eigenvectors')
% end
% 
% %==========================================================================
% % Eigenfunctions
% %==========================================================================
% [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, v, true);
% if sSimParams.b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, ...
%         sSimParams.PlotEigenFuncsM-1, mPhiAnalytic, vLambdaAnalytic, 'Analytic', [], G, ...
%         '${\bf W}_G$ kernel Eigenfunctions')
% end

%==========================================================================
% Graph-signal
%==========================================================================
if strcmp(graphSignalModel, 'Phi_c')
    [f,f_hat,coeffs] = GenerateGraphSignal(G, Phi, graphSignalModel);
elseif strcmp(graphSignalModel, 'V_c')
    [f,f_hat,coeffs] = GenerateGraphSignal(G, V, graphSignalModel);
elseif strcmp(graphSignalModel, 'alpha_K')
    [f,f_hat,coeffs] = GenerateGraphSignal(G, [], graphSignalModel);
end
%% G_tilde
%==========================================================================
% Transform
%==========================================================================
% R = eye(G.N);
% r = randperm(G.N);
% R = R(r,:);

% R = RandOrthMat(G.N);

R = 1/sqrt(G.N)*randn(G.N,G.N);

v_tilde = R*v;

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, dist_tilde, sDatasetTilde, sKernelTildeParams] = GenerateGraphTilde(v_tilde, nComponents,omega,...
    nEigs,b_normlizedLaplacian);

%==========================================================================
% Generate G_tilde
%==========================================================================
PlotGraphToGraphTransform(G, G_tilde);

%==========================================================================
% Eigenfunctions
%==========================================================================
[mPhiTildeAnalytic, vLambdaTildeAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelTildeParams, ...
    v_tilde, true);
if sSimParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, ...
        mPhiTildeAnalytic, vLambdaTildeAnalytic, 'Analytic', [], G_tilde, ...
        ['${\bf W}_{\tilde{G}}$ kernel Eigenfunctions' newline '$\tilde{v}={\bf R} v$'])
    PlotInnerProductMatrix(sSimParams, sKernelTildeParams.sDistParams, graphName, [], mPhiTildeAnalytic, 'Analytic');
end
%==========================================================================
% Eigenvectors
%==========================================================================
[V_tilde, vLambdaNumericTilde] = CalcNumericEigenvectors(nEigs, sKernelTildeParams, v_tilde);
if sSimParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, ...
        sSimParams.PlotEigenFuncsM-1, V_tilde, vLambdaNumericTilde, 'Numeric', [], G_tilde, ...
        '${\bf W}_{\tilde{G}}$ Eigenvectors')
    PlotInnerProductMatrix(sSimParams, sKernelTildeParams.sDistParams, graphName, [], V_tilde, 'Numeric');
end


%% Functional maps
% fmap: C = L2.evecs'*L2.A'*P'*L1.evecs;

if strcmp(G_tildeBasis, 'Analytic')
    [n, dim] = size(v_tilde);
    mPrTilde = diag(sKernelTildeParams.sDistParams.vPr);
    if strcmp(graphSignalModel, 'Phi_c')
        C = n^dim*mPhiTildeAnalytic'*R*mPhiTildeAnalytic;
    elseif strcmp(graphSignalModel, 'V_c')
        C = n^dim*mPhiTildeAnalytic'*mPrTilde'*R*V;
    elseif strcmp(graphSignalModel, 'alpha_K')
%         C = mPhiTildeAnalytic'*R*G.W; 
        C = n^dim*mPhiTildeAnalytic'*mPrTilde'*R*G.W;
%         C = mPhiTildeAnalytic'*mPrTilde'*R*G.W;
    end
    
elseif strcmp(G_tildeBasis, 'Numeric')
    if strcmp(graphSignalModel, 'V_c')
        C = V_tilde'*R*V;
    elseif strcmp(graphSignalModel, 'alpha_K')
        C = V_tilde'*R*G.W;
    end
else
    error('invalid basis')
end
figure; 
subplot(1,2,1)
    imagesc(C); colorbar;
    title('$C$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(1,2,2)
    imagesc(C\C); colorbar;
    title('$C^T C$', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf,'Position', [400 400 1000 400])
%% Transform
coeffsTilde = C*coeffs;

% f_str = '\sum_{i=1}^N \alpha_i K(v,v_i)';
f_str = '';

if strcmp(G_tildeBasis, 'Analytic')
    f_tilde = mPhiTildeAnalytic*coeffsTilde;
    coeffsRec = (1/n^dim)*C'*coeffsTilde;
elseif strcmp(G_tildeBasis, 'Numeric')
    f_tilde = V_tilde*coeffsTilde;
%     coeffsRec = C'*coeffsTilde;
    coeffsRec = C\coeffsTilde;
else
    error('invalid basis')
end

if strcmp(graphSignalModel, 'V_c')
    f_rec = V*coeffsRec;
elseif strcmp(graphSignalModel, 'alpha_K')
    f_rec = G.W*coeffsRec;
end

PlotGraphToGraphTransformAndBack(G, G_tilde, G, '$G$', '$\tilde{G}$', '$G$', f, f_tilde, f_rec, '$f(v)$', '$\tilde{f}(\tilde{v})$', '$f_{rec}(v)$')