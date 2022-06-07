function sPreset = GetBrazilWeatherPreset()
%% Dataset parameters
sPreset.dim                = 2;
sPreset.n                  = 148; % 148 / 196
sPreset.N                  = 296;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'BrazilWeather'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.monthNames  = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Ann'};
sDatasetParams.nLabeled    = round(0.5*sPreset.n); %20;
sPreset.sDatasetParams     = sDatasetParams;
%% Number of signals
sPreset.nSignals           = numel(sDatasetParams.monthNames);
%% Gaussian kernel width
sPreset.omega              = 7; %700; % for nystrom kernel
sPreset.omegaTilde         = 7; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 3;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20; % If M is too big eigs(W,M) is not stable numerically
sPreset.MTilde             = 450;
%% Regularizations
sPreset.gamma1             = 0.1;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.1;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 10;
%% Method parameters
sPreset.b_debugUseAnalytic = false;
sPreset.b_forceCtoIdentity = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = false;
sPreset.b_runGraphSignals  = true;
sPreset.b_maskDataFitTerm  = true;
sPreset.b_compareMethods   = true;
%% 
sPreset.dataGenTechnique = 'OneDraw';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
sPreset.sDistanceParams.earthRadius = 50; %6371; % km (approximation)
%% Other methods
sPreset.cMethods = {'RoMix', 'Rep. Thm.', 'PW', 'Nystrom', 'kNN'};
sPreset.knn = 5;
sPreset.sPwParams.regularize_epsilon = 0.01;
sPreset.sPwParams.order = 100;
end
