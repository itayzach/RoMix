function sPreset = GetBrazilWeatherPreset()
%% Dataset parameters
sPreset.dim                = 2;
sPreset.n                  = 196; % 148 / 196
sPreset.N                  = 296;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'BrazilWeather'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.monthNames  = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Ann'};
sDatasetParams.nLabeled    = sPreset.n; %20;
sPreset.sDatasetParams     = sDatasetParams;
%% Number of signals
sPreset.nSignals           = numel(sDatasetParams.monthNames);
%% Gaussian kernel width
sPreset.omega              = 1000; % for nystrom kernel
sPreset.omegaTilde         = 7; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 4;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 25*sPreset.gmmNumComponents;
sPreset.MTilde             = 25*sPreset.gmmNumComponents;
%% Regularizations
sPreset.gamma1             = 0.1;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 1e-3;
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
%% 
sPreset.dataGenTechnique = 'AddPoints';
sPreset.sDistanceParams.distType = 'Haversine'; % 'Euclidean' / 'Haversine'
sPreset.sDistanceParams.earthRadius = 6371; % km (approximation)
end