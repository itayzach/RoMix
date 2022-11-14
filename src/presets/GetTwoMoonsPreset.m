function sPreset = GetTwoMoonsPreset()
%% Dataset parameters
sPreset.dim                = 2;
sPreset.n                  = 200; %10000;
sPreset.N                  = 400; %15000;
sPreset.nLabeled           = 2;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'TwoMoons'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'RandomWalk'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sPreset.sDatasetParams     = [];
%% Number of signals
sPreset.nSignals           = 1;
%% Gaussian kernel width
sPreset.omega              = 0.3; % for nystrom kernel
sPreset.omegaTilde         = 0.3; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 8;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 20;
%% Regularizations
sPreset.gamma1             = 0.1;
sPreset.gamma2             = 1;
%% Representer theorem
sPreset.gamma1Rep          = 0.1;
sPreset.gamma2Rep          = 1;
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
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'TwoDraws';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
%% Other methods
sPreset.cMethods = {'RoMix', 'Rep. Thm.', 'VSPW', 'Nystr\"{o}m', 'w-kNN'};
sPreset.knn = 5;
sPreset.sPwParams.regularize_epsilon = 0.01;
sPreset.sPwParams.order = 100;
end
