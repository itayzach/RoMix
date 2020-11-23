%% Restart
clc; clear; close all; rng('default')

%% Parameters
%==========================================================================
% Graph data
%==========================================================================
% 'Gaussian_1D' / 'Uniform_1D' / 'Gaussian_2D' / 'Uniform_2D' / 
% 'Two_moons' /
% 'bunny' / 'twodgrid' / 'sensor' / 'minnesota' /
% 'david_seonsor' / 'swiss_roll' / 'random_ring'
graphName =  'sensor'; 

%==========================================================================
% Graph & graph-signal parameters
%==========================================================================
N           = 500;
omega       = 0.3;
estDataDist = 'Gaussian';
nComponents = 1;
b_normlizedLaplacian = true;

nEigs = 30;
samplingRatio = 0.1;
b_generateGraphSignal = true;
graphSignalModel = 'bandlimited'; % 'bandlimited' / 'V_c' / 'K_alpha'

gamma_A_eigrls = 0; 0.01;
gamma_I_eigrls = 0; 0.1; 

%==========================================================================
% Transformation
%==========================================================================
verticesTransformation = 'MDS'; % 'MDS' / 'RP' / 'RandOrth' / 'none'
tilde_gamma_A_eigrls = 0; 0.01;
tilde_gamma_I_eigrls = 0; 0.1;

%==========================================================================
% Plot parameters
%==========================================================================
if samplingRatio < 0.2
    sSimParams.b_plotSamplingPointsMarkers = true;
else
    sSimParams.b_plotSamplingPointsMarkers = false;
end
sSimParams.outputFolder               = 'figs';
sSimParams.PlotEigenFuncsM            = min(nEigs, 20);
sSimParams.b_showEigenFigures         = false;
sSimParams.b_showGraphMatricesFigures = false;
sSimParams.b_showVerticesTransform    = true;
sSimParams.b_GSPBoxPlots              = true;


% graphName = 'Two_moons';
% nEigs       = 8;
% nComponents = 2;
% gamma_A_eigrls = 0; 0.01;
% gamma_I_eigrls = 0.1; 
% tilde_gamma_A_eigrls = 0; 0.01;
% tilde_gamma_I_eigrls = 0.1;

%% Generate graph
[G, dist, sDataset, sKernelParams] = GenerateGraph(graphName, nComponents, estDataDist, omega, b_normlizedLaplacian, nEigs, N);
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
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [], G)
end

%==========================================================================
% Eigenfunctions
%==========================================================================
[mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, v, true);
if sSimParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, mPhiAnalytic, vLambdaAnalytic, 'Analytic', [], G)
end

%==========================================================================
% Eigenvalues
%==========================================================================
if sSimParams.b_showEigenFigures
    PlotSpectrum(sSimParams, sDataset, [], vLambdaAnalytic, vLambdaNumeric, []);
end

%==========================================================================
% Laplacian eigenvectors (Fourier basis)
%==========================================================================
if sSimParams.b_showEigenFigures
    G = gsp_compute_fourier_basis(G);
    PlotGraphFourierFunctions(sSimParams, G, nEigs)
end
%% Signal on the graph
if ~b_generateGraphSignal % If the dataset already contains a function
    f = sDataset.sData.y;
    sample_ind = find(sDataset.sData.y);
    f_sampled = sDataset.sData.y(sample_ind);
    f_sampled_padded = sDataset.sData.y;
else
    if ~isfield(G, 'U')
        G = gsp_compute_fourier_basis(G);
    end
    assert(strcmp(graphSignalModel, 'bandlimited'), 'this scripts supports bandlimited functions only')
    k0 = round(0.01*N);
    f_hat = zeros(N,1);
    f_hat(1:k0) = 5*sort(abs(randn(k0,1)), 'descend');
    f = GenerateGraphSignal(graphSignalModel, G.U, f_hat);
    if sSimParams.b_showEigenFigures
        f_hat = G.U*f;
        PlotGraphFourierTransform(sSimParams,G,f_hat)
    end
    
    % Sample the signal
    randPerm = randperm(G.N)';
    sample_ind = sort(randPerm(1:round(samplingRatio*G.N)));
    % non_sample_ind = setdiff(1:G.N,sample_ind);
    % sample_ind = (1:samplingRatio*G.N)';

    f_sampled = f(sample_ind);
    f_sampled_padded = zeros(size(f));
    f_sampled_padded(sample_ind) = f(sample_ind);
end


%% GSP interpolation
gsp_interp_signal = gsp_interpolate(G,f_sampled,sample_ind);

%% EigRLS interpolate (without transformation)
c = eigrls(f_sampled_padded, mPhiAnalytic, vLambdaAnalytic, gamma_A_eigrls, gamma_I_eigrls, G.L);
f_int_no_transform = mPhiAnalytic*c;

% Try with eigenVECTORS
% gamma_A_eigrls = 0; 0.01;
% gamma_I_eigrls = 0; 0.1;  
% c_eigenvecs = eigrls(sampled_signal_padded, V, vLambdaNumeric, gamma_A_eigrls, gamma_I_eigrls, G.L);
% % fprintf('c =\n\t');
% % fprintf('%f\n\t', c_eigenvecs);
% % fprintf('\n');
% eigenvecs_interp_signal = V*c_eigenvecs;

%% Plot interpolation
PlotGraphSigInterp(sSimParams, sDataset.dim, G, samplingRatio, sample_ind, f, gsp_interp_signal, f_int_no_transform, graphName);

%% Transform
if strcmp(verticesTransformation, 'none')
    return
end
% "Random projection"
if strcmp(verticesTransformation, 'RP')
    newDim = sDataset.dim;
    R = randn(sDataset.dim,newDim); % Random projection matrix
    for i = 1:sDataset.dim
       R(i,:) = R(i,:)/norm(R(i,:)); % Normalization
    end
    v_tilde = v*R; % Projection
    f_tilde_sampled_padded = f_sampled_padded;
elseif strcmp(verticesTransformation, 'RandOrth')
    % % Linear combination of Gaussians
    % R = randn(G.N, G.N);
    % % for i = 1:G.N
    % %    R(:,i) = R(:,i)/norm(R(:,i)); % Normalization
    % % end
    % for i = 1:G.N
    %    R(i,:) = R(i,:)/norm(R(i,:)); % Normalization
    % end
    % v_tilde = R*v;

    % Random orthonormal matrix
    R = RandOrthMat(G.N);
    v_tilde = R*v;
    f_tilde_sampled_padded = R*f_sampled_padded;

    % % Whitening matrix
    % [v_tilde, R, invR] = whiten(v,1e-5);
    
elseif strcmp(verticesTransformation, 'MDS')
    [v_tilde, mdsLambda] = cmdscale(dist);
%     [v_tilde, mdsLambda] = mdscale(dist, sDataset.dim);
%     Rf = f;
    f_tilde_sampled_padded = f_sampled_padded;
end


%% Plot vertices transformation
if sSimParams.b_showVerticesTransform
    PlotVerticesTransformation(sSimParams,sDataset.dim,v,v_tilde,graphName,verticesTransformation);
end

%% Generate G_tilde
[G_tilde, dist_tilde, sDatasetTilde, sKernelParamsTilde] = GenerateGraphTilde(v_tilde, nComponents,omega,nEigs,b_normlizedLaplacian);

%==========================================================================
% Eigenfunctions
%==========================================================================
[mPhiTildeAnalytic, vLambdaTildeAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParamsTilde, v_tilde, true);
if sSimParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v_tilde, [], 0, sSimParams.PlotEigenFuncsM-1, mPhiTildeAnalytic, vLambdaTildeAnalytic, 'Analytic', [], G_tilde)
end

%==========================================================================
% Eigenvectors
%==========================================================================
[~, vLambdaTildeNumeric] = CalcNumericEigenvectors(nEigs, sKernelParamsTilde, v_tilde);
if sSimParams.b_showEigenFigures
    PlotSpectrum(sSimParams, sDatasetTilde, [], vLambdaTildeAnalytic, vLambdaTildeNumeric, []);
end

%% Distances
fprintf('distances diffs norm = %d\n', norm(dist_tilde - dist, 'fro'));

%% Find coeffs for graph-signal - EigRLS
c_tilde = eigrls(f_tilde_sampled_padded, mPhiTildeAnalytic, vLambdaTildeAnalytic, tilde_gamma_A_eigrls, tilde_gamma_I_eigrls, G.L);
% c_tilde = lsqminnorm(mPhiTildeAnalytic,R*f);
% c_tilde = (mPhiTildeAnalytic.'*mPhiTildeAnalytic)\(mPhiTildeAnalytic.'*(R*f));
% c_tilde = pinv(mPhiTildeAnalytic)*R*f;

f_tilde_int_eigrls = mPhiTildeAnalytic*c_tilde; % Interpolate
if strcmp(verticesTransformation, 'RandOrth')
    f_int_eigrls = R'*f_tilde_int_eigrls;
else
    f_int_eigrls = f_tilde_int_eigrls;
end

PlotGraphSigInterpTransform(sSimParams, sDataset.dim, G, f, G_tilde, f_tilde_int_eigrls, samplingRatio, sample_ind, graphName, verticesTransformation, 'EigRLS');
% PlotGraphSigInterpTransformNoReturn(sSimParams, sDataset.dim, G, f, G_tilde, f_tilde_int_eigrls, samplingRatio, sample_ind, graphName, verticesTransformation, 'EigRLS');


%% Find coeffs for graph-signal - fminsearch
% maxIter_fminsearch = 100;
% fun = @(c)f_sampled'*G.L(sample_ind,sample_ind)*f_sampled - c'*mPhiTildeAnalytic.'*G_tilde.L*mPhiTildeAnalytic*c;
% c0 = randn(nEigs,1);
% options = optimset('MaxIter',maxIter_fminsearch);
% c_tilde = fminsearch(fun,c0,options);
% f_tilde_int_fminsearch = mPhiTildeAnalytic*c_tilde;  % Interpolate
% if strcmp(verticesTransformation, 'RandOrth')
%     f_int_fminsearch = R'*f_tilde_int_fminsearch;
% else
%     f_int_fminsearch = f_tilde_int_fminsearch;
% end

% PlotGraphSigInterpTransform(sSimParams, sDataset.dim, G, f, G_tilde, f_tilde_int_fminsearch, samplingRatio, sample_ind, graphName, verticesTransformation, 'fminsearch');
% PlotGraphSigInterpTransformNoReturn(sSimParams, sDataset.dim, G, f, G_tilde, f_tilde_int_fminsearch, samplingRatio, sample_ind, graphName, verticesTransformation, 'fminsearch');
