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
graphSignalModel = 'V_c'; % 'alpha_K' / 'V_c'

nEigs = 30;
samplingRatio = 0;
nysRatio = 0.8;
%==========================================================================
% Transformation
%==========================================================================
verticesTransformation = 'MDS'; % 'MDS' / 'RP' / 'RandOrth'

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
sSimParams.b_showEigenFigures         = true;
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
% [f,f_hat,alpha] = GenerateGraphSignal(G, [], graphSignalModel);
[f,f_hat,a] = GenerateGraphSignal(G, V, graphSignalModel);
%% G_tilde
%==========================================================================
% Transform
%==========================================================================
R = eye(G.N);
r = randperm(G.N);
R = R(r,:);
% R = RandOrthMat(G.N);
v_tilde = R*v;

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, dist_tilde, sDatasetTilde, sKernelParamsTilde] = GenerateGraphTilde(v_tilde, nComponents,omega,...
    nEigs,b_normlizedLaplacian);

%==========================================================================
% Generate G_tilde
%==========================================================================
PlotGraphToGraphTransform(G, G_tilde);
% %==========================================================================
% % Eigenfunctions
% %==========================================================================
% [mPhiTildeAnalytic, vLambdaTildeAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParamsTilde, ...
%     v_tilde, true);
% if sSimParams.b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, ...
%         mPhiTildeAnalytic, vLambdaTildeAnalytic, 'Analytic', [], G_tilde, ...
%         ['${\bf W}_{\tilde{G}}$ kernel Eigenfunctions' newline '$\tilde{v}={\bf R} v$'])
% end

%==========================================================================
% Eigenvectors
%==========================================================================
[V_tilde, vLambdaNumericTilde] = CalcNumericEigenvectors(nEigs, sKernelParamsTilde, v_tilde);
if sSimParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, ...
        sSimParams.PlotEigenFuncsM-1, V_tilde, vLambdaNumericTilde, 'Numeric', [], G_tilde, ...
        '${\bf W}_{\tilde{G}}$ Eigenvectors')
end
%% fmap attempt
% fmap: C = L2.evecs'*L2.A'*P'*L1.evecs;
% C = mPhiTildeAnalytic'*R*G.W;
% C = V_tilde'*R*G.W;
C = V_tilde'*R*V;
figure; imagesc(C); colorbar;
%% Transform
% mPhiTransformed = R*mPhiTildeAnalytic;
% if sSimParams.b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, ...
%         mPhiTransformed, vLambdaNumeric, 'Analytic', [], G, ...
%         ['${\bf W}_G$ kernel Eigenfunctions' newline 'after interpolation'])
% end
% PlotGraphToGraphTransform(G, G_tilde, mPhiTransformed(:,1), mPhiTildeAnalytic(:,1), '$\phi_0(v) = {\bf R} \tilde{\phi}_0(\tilde{v})$', '$\tilde{\phi}_0(\tilde{v})$');
% 
% 
% mPhiTransformed = R'*mPhiTildeAnalytic;
% if sSimParams.b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, ...
%         mPhiTransformed, vLambdaNumeric, 'Analytic', [], G, ...
%         ['${\bf W}_G$ kernel Eigenfunctions' newline 'after interpolation'])
% end
% PlotGraphToGraphTransform(G, G_tilde, mPhiTransformed(:,1), mPhiTildeAnalytic(:,1), '$\phi_0(v) = {\bf R}^T \tilde{\phi}_0(\tilde{v})$', '$\tilde{\phi}_0(\tilde{v})$');
% 
% mPhiTransformed = mPhiTildeAnalytic;
% if sSimParams.b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, ...
%         mPhiTransformed, vLambdaNumeric, 'Analytic', [], G, ...
%         ['${\bf W}_G$ kernel Eigenfunctions' newline 'after interpolation'])
% end
% PlotGraphToGraphTransform(G, G_tilde, mPhiTransformed(:,1), mPhiTildeAnalytic(:,1), '$\phi_0(v) = \tilde{\phi}_0(\tilde{v})$', '$\tilde{\phi}_0(\tilde{v})$');
b = C*a;
f_tilde = V_tilde*b;
% f_str = '\sum_{i=1}^N \alpha_i K(v,v_i)';
f_str = '';
PlotGraphToGraphTransform(G, G_tilde, f, f_tilde, ['$f(v) = ' f_str '$'], '$\tilde{f}(\tilde{v}) = ??? $');


