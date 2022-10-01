function sPreset = GetBulgariBeaconsPreset()
%% Dataset parameters
sPreset.dim                = 9;
sPreset.NLat               = 100; % 100 % 120
sPreset.NLon               = 100; % 100 % 120
sPreset.N                  = sPreset.NLat*sPreset.NLon;
sPreset.n                  = 6000;
sPreset.nLabeled           = [25 50 100 150];
sPreset.k                  = round(0.1*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'BulgariBeacons'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.xTickNames  = cellfun(@num2str,(num2cell(sPreset.nLabeled)), 'UniformOutput', false);
sDatasetParams.b_zoomInAcc = true;
sPreset.sDatasetParams     = sDatasetParams;
%% Number of signals
sPreset.nSignals           = 1;
%% Gaussian kernel width
sPreset.omega              = 1.5; % for nystrom kernel
sPreset.omegaTilde         = 1.5; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 100;
sPreset.gmmNumComponents   = 8;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 10; % If M is too big eigs(W,M) is not stable numerically
sPreset.MTilde             = 100;
%% Regularizations
sPreset.gamma1             = 0.05;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.05;
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
%% Other methods
sPreset.cMethods = {'RoMix', 'Rep. Thm.', 'PW', 'Nystrom', 'kNN'};
sPreset.knn = 5;
sPreset.sPwParams.regularize_epsilon = 0;
sPreset.sPwParams.order = 50;
end
