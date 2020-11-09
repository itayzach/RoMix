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
G_basisModel     = 'Numeric';  % 'Analytic' / 'Numeric' / 'LaplacianEigenvectors' / 'representerBasis'
Gtilde_basisModel = 'Analytic'; % 'Analytic' / 'Numeric' / 'LaplacianEigenvectors'

verticesTransform = 'randomMatrix'; % 'permutation' / 'RandOrthMat' / 'randomMatrix' / 'eye'
funcTransform     = 'pinv(Btilde)RB'; % 'pinv(Btilde)RB' / 'pinv(Btilde)B'
G_nEigs = 30;
Gtilde_nEigs = 30;
samplingRatio = 1;0.8;


% % TODO: the following does not work:
% graphSignalModel = 'U_fhat';
% G_basisModel     = 'LaplacianEigenvectors';
% Gtilde_basisModel = 'Numeric';
% verticesTransform = 'randomMatrix';
% funcTransform     = 'pinv(Btilde)RB'; 
%==========================================================================
% Plot parameters
%==========================================================================
if samplingRatio < 0.2
    sSimParams.b_plotSamplingPointsMarkers = true;
else
    sSimParams.b_plotSamplingPointsMarkers = false;
end
% sSimParams.b_plotSamplingPointsMarkers = false;

sSimParams.outputFolder               = 'figs';
sSimParams.PlotEigenFuncsM            = min(G_nEigs, 12);
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = false;
sSimParams.b_GSPBoxPlots              = false;
sSimParams.b_plotTransformation       = false;
sSimParams.b_calcCoeffsInterpolated   = false;
sSimParams.b_interpolateOnGWithLS     = true;
%% Generate graph
[G, dist, sDataset, sKernelParams] = GenerateGraph(graphName, nComponents, ...
    estDataDist, omega, b_normlizedLaplacian, G_nEigs);
v_N = sDataset.sData.x; % For convenience

if sSimParams.b_showGraphMatricesFigures
    PlotGraphMatrices(G, b_normlizedLaplacian);
end
if sDataset.dim == 1
    sSimParams.b_GSPBoxPlots = false;
end

%%
%==========================================================================
% Sampling matrix
%==========================================================================
randPerm = randperm(G.N)';
sampleInd = sort(randPerm(1:round(samplingRatio*G.N)));
S = zeros(G.N,1);
S(sampleInd) = 1;
S = diag(S);

v_n = v_N(sampleInd,:);
n = length(v_n);
%% Basis on G
%==========================================================================
% Eigenfunctions
%==========================================================================
if strcmp(G_basisModel, 'Analytic')
    [mPhi, vLambdaAnalytic] = CalcAnalyticEigenfunctions(G_nEigs, sKernelParams, v_n, true);
    if sSimParams.b_showEigenFigures
        PlotEigenfuncvecScatter(sSimParams, estDataDist, v_n, [], 0, ...
            sSimParams.PlotEigenFuncsM-1, mPhi, vLambdaAnalytic, 'Analytic', [], G, ...
            '${\bf W}_G$ kernel Eigenfunctions')
    end
    Bn = mPhi;
    Bn_eigvals = vLambdaAnalytic;
%==========================================================================
% Eigenvectors
%==========================================================================    
elseif strcmp(G_basisModel, 'Numeric')
    [V, vLambdaNumeric] = CalcNumericEigenvectors(G_nEigs, sKernelParams, v_n);
    if sSimParams.b_showEigenFigures
        PlotEigenfuncvecScatter(sSimParams, estDataDist, v_n, [], 0, ...
            sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [], G, ...
            '${\bf W}_G$ Eigenvectors')
    end
    Bn = V;
    Bn_eigvals = vLambdaNumeric;
    
    [V, vLambdaNumeric] = CalcNumericEigenvectors(G_nEigs, sKernelParams, v_N);
    if sSimParams.b_showEigenFigures
        PlotEigenfuncvecScatter(sSimParams, estDataDist, v_N, [], 0, ...
            sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [], G, ...
            '${\bf W}_G$ Eigenvectors')
    end
    BN = V;
    BN_eigvals = vLambdaNumeric;
%==========================================================================
% Laplacian Eigenvectors
%==========================================================================    
elseif strcmp(G_basisModel, 'LaplacianEigenvectors')
    G = gsp_compute_fourier_basis(G); 
    Bn = G.U;
    Bn_eigvals = G.e;
else
    error('invalid basis')
end



%%
%==========================================================================
% Graph-signal
%==========================================================================
if strcmp(graphSignalModel, 'Phi_c')
    f_title = '$f(v) = {\bf \Phi} c$';
elseif strcmp(graphSignalModel, 'V_c')
    f_title = '$f(v) = {\bf V} \alpha$';
elseif strcmp(graphSignalModel, 'K_alpha')
    f_title = '$f(v) = {\bf K} \alpha$';
elseif strcmp(graphSignalModel, 'U_fhat')
    f_title = '$f(v) = {\bf U} \hat{f}$';
end
[f, coeffsGroundTruthForDebug] = GenerateGraphSignal(Bn, graphSignalModel);

% f_sampled_padded = zeros(size(f));
% f_sampled_padded(sampleInd) = f(sampleInd);

%==========================================================================
% Coefficients
%==========================================================================
if sSimParams.b_interpolateOnGWithLS
    % find coeffs by projecting f onto V ( c = V'*f ) and use least squares since f is sampled
    coeffsLS = pinv(Bn)*f;
%     error('TODO: f_interp_no_transform is not interpolated! it is 800x1 instead of 1000x1')
    f_interp_no_transform = Bn*coeffsLS;
    f_interp_no_transform_title = '$f_{{\bf int}}(v) = {\bf B}({\bf SB})^\dagger f_{\mathcal{S}}(v)$';
    PlotGraphToGraphTransform(sSimParams, G, '$G$',  G, '$G$', [], [], ...
        f, f_title, f_interp_no_transform, f_interp_no_transform_title, [], [], sampleInd)
end
% norm(f - f_interp_no_transform)/norm(f)
%% G_tilde
%==========================================================================
% Transform
%==========================================================================
if strcmp(verticesTransform, 'permutation')
    R = eye(G.N);
    r = randperm(G.N);
    R = R(r,:);
elseif strcmp(verticesTransform, 'eye')
    R = eye(G.N);
elseif strcmp(verticesTransform, 'RandOrthMat')
    R = RandOrthMat(G.N);
elseif strcmp(verticesTransform, 'randomMatrix')
    R = (1/sqrt(G.N))*randn(G.N,G.N);
end
v_tilde = R*v_N; % TODO: should I use all v or just v_n?

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, dist_tilde, sDatasetTilde, sKernelTildeParams] = GenerateGraphTilde(v_tilde, nComponents,omega,...
    Gtilde_nEigs,b_normlizedLaplacian);

%==========================================================================
% Eigenfunctions
%==========================================================================
if strcmp(Gtilde_basisModel, 'Analytic')
    [mPhiTilde, vLambdaAnalyticTilde] = CalcAnalyticEigenfunctions(Gtilde_nEigs, sKernelTildeParams, ...
        v_tilde, true);
    if sSimParams.b_showEigenFigures
        PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, ...
            mPhiTilde, vLambdaAnalyticTilde, 'Analytic', [], G_tilde, ...
            ['${\bf W}_{\tilde{G}}$ kernel Eigenfunctions' newline '$\tilde{v}={\bf R} v$'])
        PlotInnerProductMatrix(sSimParams, sKernelTildeParams.sDistParams, graphName, [], mPhiTilde, 'Analytic');
    end
    Btilde = mPhiTilde;
    Btilde_eigvals = vLambdaAnalyticTilde;
    
elseif strcmp(Gtilde_basisModel, 'Numeric')
%==========================================================================
% Eigenvectors
%==========================================================================
    [V_tilde, vLambdaNumericTilde] = CalcNumericEigenvectors(Gtilde_nEigs, sKernelTildeParams, v_tilde);
    if sSimParams.b_showEigenFigures
        PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, ...
            sSimParams.PlotEigenFuncsM-1, V_tilde, vLambdaNumericTilde, 'Numeric', [], G_tilde, ...
            '${\bf W}_{\tilde{G}}$ Eigenvectors')
        PlotInnerProductMatrix(sSimParams, sKernelTildeParams.sDistParams, graphName, [], V_tilde, 'Numeric');
    end
    Btilde = V_tilde;
    Btilde_eigvals = vLambdaNumericTilde;
%==========================================================================
% Laplacian Eigenvectors
%==========================================================================    
elseif strcmp(Gtilde_basisModel, 'LaplacianEigenvectors')
    G_tilde = gsp_compute_fourier_basis(G_tilde); 
    Btilde = G_tilde.U;
    Btilde_eigvals = G_tilde.e;
end
%% Functional maps: C = L2.evecs'*L2.A'*P'*L1.evecs;
%TODO: try with mPrTidle? C = mPhiTilde'*mPrTilde'*R*V;
if strcmp(Gtilde_basisModel, 'Analytic')
    [n, dim] = size(v_tilde);
    mPrTilde = diag(sKernelTildeParams.sDistParams.vPr);
end

if strcmp(funcTransform, 'pinv(Btilde)RB')
    C = pinv(Btilde)*R*Bn;
elseif strcmp(funcTransform, 'pinv(Btilde)B')
    C = pinv(Btilde)*Bn;
else
    error('select funcTransform');
end
T = Btilde*C*pinv(Bn);
TS = Btilde*C*pinv(S*Bn);

if sSimParams.b_plotTransformation
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
end
%% Transform
% coeffsTilde = C*coeffsLS;

% TODO: Make following lines are be like S*T*f_sampled_padded;
f_tilde_sampled = TS(sampleInd,sampleInd)*f;
f_tilde_sampled_padded = zeros(size(f));
f_tilde_sampled_padded(sampleInd) = f_tilde_sampled;


if strcmp(Gtilde_basisModel, 'Analytic')
    f_tilde_interp_title = '$\tilde{f}_{{\bf int}}(\tilde{v}) = \tilde{{\bf \Phi}} \tilde{c}$';
elseif strcmp(Gtilde_basisModel, 'Numeric')
    f_tilde_interp_title = '$\tilde{f}_{{\bf int}}(\tilde{v}) = \tilde{{\bf V}} \tilde{\alpha}$';
elseif strcmp(Gtilde_basisModel, 'LaplacianEigenvectors')
    f_tilde_interp_title = '$\tilde{f}_{{\bf int}}(\tilde{v}) = \tilde{{\bf U}} \tilde{\hat{f}}$';
else
    error('invalid basis')
end

% coeffsTilde = ( (S*Btilde).' * (S*Btilde) ) \ ((S*Btilde).' * f_tilde_sampled_padded);
coeffsTilde = eigrls(f_tilde_sampled_padded, Btilde, Btilde_eigvals, 0, 0, G_tilde.L);
f_tilde_interp = Btilde*coeffsTilde;

if sSimParams.b_calcCoeffsInterpolated
    % TODO: multiply by (1/samplingRatio)?
    % TODO: can we use C' instead of inv(C)?
    coeffs_interp = C\coeffsTilde; 
    figure; 
    scatter(1:length(coeffsLS), coeffsLS, 'bx', 'DisplayName', 'Least-squares')
    hold on
    scatter(1:length(coeffs_interp), coeffs_interp, 'ro', 'DisplayName', 'Interpolated')
    legend()
end

f_interp = pinv(T)*f_tilde_interp;
f_interp_title = '$f_{{\bf int}}(v) = {\bf T}^\dagger \tilde{f}_{{\bf int}}(\tilde{v})$';

PlotGraphToGraphTransform(sSimParams, G, '$G$', G_tilde, '$\tilde{G}$', ...
    G, '$G$', f, f_title, f_tilde_interp, f_tilde_interp_title, f_interp, f_interp_title, sampleInd)


% norm(f - f_interp)/norm(f)