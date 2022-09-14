function sPreset = GetSwissRollPreset(b_randn, b_gmmLatent)
if ~exist('b_randn', 'var')
    b_randn = false;
end
if ~exist('b_gmmLatent','var')
    b_gmmLatent = false;
end
%% Dataset parameters
sPreset.dim                = 3*(~b_gmmLatent) + 2*b_gmmLatent;
sPreset.n                  = 2000;
sPreset.nLabeled           = 1500;
sPreset.N                  = 5000;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'SwissRoll'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'NormLap'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.a           = 1;
sDatasetParams.maxTheta    = 4*pi;
sDatasetParams.height      = 20;
sDatasetParams.b_randn     = b_randn;
sDatasetParams.b_gmmLatent = b_gmmLatent;
if b_gmmLatent
    gmmNumComponents = 5;
else
    gmmNumComponents = 20;
end
if b_randn
    maxS = SwissRollArclength(sDatasetParams.maxTheta);
    for c = 1:gmmNumComponents
        sDatasetParams.sigma{c}    = [maxS/(5*gmmNumComponents), sDatasetParams.height/5];
        sDatasetParams.mu{c}       = [(c/(gmmNumComponents+1))*maxS, sDatasetParams.height/2];
        sDatasetParams.compProp{c} = 1/gmmNumComponents;
    end
end
sPreset.sDatasetParams     = sDatasetParams;
%% Number of signals
sPreset.nSignals           = 1;
%% Gaussian kernel width
% eps value was taken from:
% "Parsimonious representation of nonlinear dynamical systems through manifold learning: 
% A chemotaxis case study"
% where K(x_i,x_j) = exp(||x_i-x_j||^2/eps^2)
eps                        = sqrt(5);
sPreset.omega              = eps/sqrt(2); % for nystrom kernel
sPreset.omegaTilde         = eps/sqrt(2); % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-5;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = gmmNumComponents;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 50*sPreset.gmmNumComponents*(~b_randn) + sPreset.M*b_randn;
%% Regularizations
sPreset.gamma1             = 1e-5;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 1e-5;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 10;
%% Method parameters
sPreset.b_debugUseAnalytic = true;
sPreset.b_forceCtoIdentity = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = true;
sPreset.b_runGraphSignals  = true;
sPreset.b_maskDataFitTerm  = false;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'OneDraw';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end
