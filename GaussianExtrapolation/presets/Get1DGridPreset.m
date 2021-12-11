function sPreset = Get1DGridPreset()
%% Dataset parameters
sPreset.dim                = 1;
sPreset.n                  = 1500;
sPreset.N                  = 3000;
sPreset.k                  = 3;
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'Grid'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'RandomWalk'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.xMin        = -1*ones(sPreset.dim,1);
sDatasetParams.xMax        = 1*ones(sPreset.dim,1);
sPreset.sDatasetParams     = sDatasetParams;
%% Gaussian kernel width
L = sDatasetParams.xMax(1) - sDatasetParams.xMin(1);
sPreset.omega              = 2*sqrt(L/sPreset.n); % for nystrom kernel
sPreset.omegaTilde         = 2*sqrt(L/sPreset.n); % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 1;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 30;
sPreset.gamma1             = 0.001;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.001;
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