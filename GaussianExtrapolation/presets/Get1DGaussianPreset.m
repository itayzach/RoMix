function sPreset = Get1DGaussianPreset()
%% Dataset parameters
sPreset.dim                = 1;
sPreset.n                  = 1500;
sPreset.N                  = 3000;
sPreset.k                  = 3;
sPreset.nGenDataCompnts    = 1;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'Gaussian'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
for c = 1:sPreset.nGenDataCompnts
    sDatasetParams.mu{c}    = 10*(c-1)*ones(1,sPreset.dim);
    sDatasetParams.sigma{c} = 1*eye(sPreset.dim);
end
sPreset.sDatasetParams     = sDatasetParams;
%% Number of signals
sPreset.nSignals           = 1;
%% Gaussian kernel width
sPreset.omega              = 0.3; % for nystrom kernel
sPreset.omegaTilde         = 0.3; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = sPreset.nGenDataCompnts;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 30;
%% Regularizations
sPreset.gamma1             = 0.001;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.001;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 1;
%% Method parameters
sPreset.b_debugUseAnalytic = false;
sPreset.b_forceCtoIdentity = false;
sPreset.b_normalizePhi     = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = false;
sPreset.b_runGraphSignals  = true;
sPreset.b_maskDataFitTerm  = false;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'AddPoints';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end
