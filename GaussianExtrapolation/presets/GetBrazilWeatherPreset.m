function sPreset = GetBrazilWeatherPreset()
%% Dataset parameters
sPreset.dim                = 2;
sPreset.n                  = 246;
sPreset.N                  = 296;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'BrazilWeather'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sPreset.sDatasetParams     = [];
%% Gaussian kernel width
sPreset.omega              = 7; % for nystrom kernel
sPreset.omegaTilde         = 7; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 4;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 25*sPreset.gmmNumComponents;
sPreset.MTilde             = 25*sPreset.gmmNumComponents;
sPreset.gamma1             = 0.1;
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
sPreset.b_interpEigenvecs  = false;
sPreset.b_runGraphSignals  = true;
sPreset.b_maskDataFitTerm  = false;
%% 
sPreset.dataGenTechnique = 'AddPoints';
sPreset.distType = 'Haversine'; % 'Euclidean' / 'Haversine'

end