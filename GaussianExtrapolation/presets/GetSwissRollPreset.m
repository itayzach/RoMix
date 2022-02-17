function sPreset = GetSwissRollPreset()
%% Dataset parameters
sPreset.dim                = 3;
sPreset.n                  = 4000;
sPreset.N                  = 5000;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'SwissRoll'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'RandomWalk'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sPreset.sDatasetParams     = [];
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
sPreset.gmmNumComponents   = 20;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 10;
sPreset.MTilde             = 50*sPreset.gmmNumComponents;
%% Regularizations
sPreset.gamma1             = 1e-5;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 1e-5;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 10;
%% Method parameters
sPreset.b_debugUseAnalytic = false;
sPreset.b_forceCtoIdentity = false;
sPreset.b_normalizePhi     = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = true;
sPreset.b_runGraphSignals  = true;
sPreset.b_maskDataFitTerm  = false;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'AddPoints';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end
