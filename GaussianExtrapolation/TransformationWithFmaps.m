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
graphSignalModel = 'V_c'; % 'K_alpha' / 'V_c' / 'Phi_c' / 'U_fhat'
G_tildeBasis = 'Analytic'; % 'Analytic' / 'Numeric'
verticesTransform = 'randomMatrix'; % 'permutation' / 'RandOrthMat' / 'randomMatrix'
nEigs = 30;
samplingRatio = 0.05;
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
    [f,f_hat,coeffsGroundTruthForDebug] = GenerateGraphSignal(G, Phi, graphSignalModel);
elseif strcmp(graphSignalModel, 'V_c')
    [f,f_hat,coeffsGroundTruthForDebug] = GenerateGraphSignal(G, V, graphSignalModel);
elseif strcmp(graphSignalModel, 'K_alpha')
    [f,f_hat,coeffsGroundTruthForDebug] = GenerateGraphSignal(G, [], graphSignalModel);
elseif strcmp(graphSignalModel, 'U_fhat')
    if ~isfield(G, 'U')
        G = gsp_compute_fourier_basis(G);
    end
    [f,f_hat,coeffsGroundTruthForDebug] = GenerateGraphSignal(G, [], graphSignalModel);
end
randPerm = randperm(G.N)';
sampleInd = sort(randPerm(1:round(samplingRatio*G.N)));
f_sampled = f(sampleInd);
f_sampled_padded = zeros(size(f));
f_sampled_padded(sampleInd) = f(sampleInd);

%==========================================================================
% Coefficients
%==========================================================================
% find coeffs by projecting f onto V ( c = V'*f ) and use least squares since f is sampled
J = diag(sign(abs(f_sampled_padded)));
if strcmp(graphSignalModel, 'Phi_c')
    coeffsLS = pinv(J*mPhiAnalytic)*f_sampled_padded;
elseif strcmp(graphSignalModel, 'V_c')
    coeffsLS = pinv(J*V)*f_sampled_padded;
elseif strcmp(graphSignalModel, 'K_alpha')
    warning('is K a basis? can I do that?');
    coeffsLS = pinv(J*G.W)*f_sampled_padded;
elseif strcmp(graphSignalModel, 'U_fhat')
    coeffsLS = pinv(J*G.U)*f_sampled_padded; % same as iGFT: coeffs = f_hat = G.U'*f
end


%% G_tilde
%==========================================================================
% Transform
%==========================================================================
if strcmp(verticesTransform, 'permutation')
    R = eye(G.N);
    r = randperm(G.N);
    R = R(r,:);
elseif strcmp(verticesTransform, 'RandOrthMat')
    R = RandOrthMat(G.N);
elseif strcmp(verticesTransform, 'randomMatrix')
    R = (1/sqrt(G.N))*randn(G.N,G.N);
end
v_tilde = R*v;

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, dist_tilde, sDatasetTilde, sKernelTildeParams] = GenerateGraphTilde(v_tilde, nComponents,omega,...
    nEigs,b_normlizedLaplacian);

%==========================================================================
% Generate G_tilde
%==========================================================================
if size(v,2) > 1
%     PlotGraphToGraphTransform(G, G_tilde);
    PlotGraphToGraphTransform(sSimParams, G, '$G$', G_tilde, '$\tilde{G}$');

end

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
        C = mPhiTildeAnalytic'*R*mPhiTildeAnalytic;
    elseif strcmp(graphSignalModel, 'V_c')
%         C = mPhiTildeAnalytic'*mPrTilde'*R*V;
        C = mPhiTildeAnalytic'*R*V;
    elseif strcmp(graphSignalModel, 'K_alpha')
%         C = mPhiTildeAnalytic'*R*G.W; 
%         C = n^dim*mPhiTildeAnalytic'*mPrTilde'*R*G.W;
%         C = mPhiTildeAnalytic'*mPrTilde'*R*G.W;
%         C = mPhiTildeAnalytic'*mPrTilde'*R*G.W;
        C = mPhiTildeAnalytic'*R*G.W;
    elseif strcmp(graphSignalModel, 'U_fhat')
        C = mPhiTildeAnalytic'*R*G.U;
    end
    
elseif strcmp(G_tildeBasis, 'Numeric')
    if strcmp(graphSignalModel, 'V_c')
        C = V_tilde'*R*V;
    elseif strcmp(graphSignalModel, 'K_alpha')
        C = V_tilde'*R*G.W;
    elseif strcmp(graphSignalModel, 'U_fhat')
        C = V_tilde'*R*G.U;        
    end
else
    error('invalid basis')
end
figure; 
subplot(1,3,1)
    imagesc(C); colorbar;
    title('$C$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(1,3,2)
    imagesc(C\C); colorbar;
    title('$C^{-1} C$', 'Interpreter', 'latex', 'FontSize', 14)
subplot(1,3,3)
    imagesc(C'*C); colorbar;
    title('$C^T C$', 'Interpreter', 'latex', 'FontSize', 14)    
set(gcf,'Position', [400 400 1500 400])
%% Transform
coeffsTilde = C*coeffsLS; 
coeffsTildeDebug = C*coeffsGroundTruthForDebug;


if strcmp(G_tildeBasis, 'Analytic')
    f_tilde = mPhiTildeAnalytic*coeffsTilde;
    coeffsRec = C\coeffsTilde;
    f_tilde_title = '$\tilde{f}(\tilde{v}) = \tilde{{\bf \Phi}} \tilde{c}$';
elseif strcmp(G_tildeBasis, 'Numeric')
    f_tilde = V_tilde*coeffsTilde;
%     coeffsRec = C'*coeffsTilde;
    coeffsRec = C\coeffsTilde;
    f_tilde_title = '$\tilde{f}(\tilde{v}) = \tilde{{\bf V}} \tilde{c}$';
else
    error('invalid basis')
end

if strcmp(graphSignalModel, 'V_c')
    f_interp = V*coeffsRec;
    f_interp_title = '$f_{{\bf int}}(v) = {\bf V} c_{{\bf int}}$';
    f_title = '$f(v) = {\bf V} c$';
elseif strcmp(graphSignalModel, 'K_alpha')
    f_interp = G.W*coeffsRec;
    f_interp_title = '$f_{{\bf int}}(v) = {\bf K} \alpha_{{\bf int}}$';
    f_title = '$f(v) = {\bf K} \alpha$';
elseif strcmp(graphSignalModel, 'U_fhat')
    f_interp = G.U*coeffsRec;
    f_interp_title = '$f_{{\bf int}}(v) = {\bf U} \hat{f}_{{\bf int}}$';
    f_title = '$f(v) = {\bf U} \hat{f}$';    
end

figure; 
scatter(1:length(coeffsGroundTruthForDebug), coeffsGroundTruthForDebug, 'ro')
hold on
scatter(1:length(coeffsLS), coeffsLS, 'bx')
scatter(1:length(coeffsRec), coeffsRec, 'k+')
legend('Ground-truth', 'Least-squares', 'Interpolated')

PlotGraphToGraphTransform(sSimParams, G, '$G$', G_tilde, '$\tilde{G}$', G, '$G$', f, f_title, f_tilde, f_tilde_title, f_interp, f_interp_title, sampleInd)
