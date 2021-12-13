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
    sDatasetParams.mu{c}    = 3*(c-1)*ones(1,sPreset.dim);
    sDatasetParams.sigma{c} = 1*eye(sPreset.dim);
end
sPreset.sDatasetParams     = sDatasetParams;
%% Gaussian kernel width
sPreset.omega              = 0.3; % for nystrom kernel
sPreset.omegaTilde         = 0.3; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 1;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 30;
%% Regularizations
sPreset.gamma1             = 0.001;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.001;
%% Number of runs (=realizations)
sPreset.R                  = 1;
%% Method parameters
sPreset.b_debugUseAnalytic = false;
sPreset.b_forceCtoIdentity = false;
sPreset.b_normalizePhi     = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = true;
sPreset.b_runGraphSignals  = false;
sPreset.b_maskDataFitTerm  = false;
%% 
sPreset.dataGenTechnique = 'AddPoints';
sPreset.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end