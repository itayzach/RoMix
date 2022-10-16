function sPreset = GetCoraVgaePreset()
%% VAE parameters
b_runVAE = true;
latentDim                  = 16;
epochs                     = 200;
batchSize                  = 128;
b_forceLoadTrainedVAE      = true; % Forces load of a trained model on entire MNIST (60,000 train)
%% Dataset parameters
sPreset.dim                = latentDim;
sPreset.n                  = 2700; % cannot be less than 28*28 due to GMM
sPreset.N                  = 2708;
sPreset.nLabeled           = round(0.1*sPreset.n);
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'CoraLatentVGAE'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.b_runVAE    = b_runVAE;
sDatasetParams.latentDim   = latentDim;
sDatasetParams.epochs      = epochs;
sDatasetParams.batchSize   = batchSize;
sDatasetParams.b_forceLoadTrainedVAE = b_forceLoadTrainedVAE;
sPreset.sDatasetParams     = sDatasetParams;
assert(sPreset.nLabeled <= sPreset.n)
%% Number of signals
sPreset.nSignals           = 1; % After conversion from number of classes
%% Gaussian kernel width
sPreset.omega              = 0.1*sPreset.dim; % for nystrom kernel
sPreset.omegaTilde         = 0.1*sPreset.dim; % for our method
%% GMM params
sPreset.gmmRegVal          = 0.1;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 1;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 500;
%% Regularizations
sPreset.gamma1             = 0.1;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.1;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 1;
%% Method parameters
sPreset.b_debugUseAnalytic = false;
sPreset.b_forceCtoIdentity = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = true;
sPreset.b_runGraphSignals  = false;
sPreset.b_maskDataFitTerm  = true;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'OneDraw';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
%% Other methods
sPreset.cMethods = {'RoMix', 'Rep. Thm.', 'VSPW', 'Nystr\"{o}m', 'w-kNN'};
sPreset.knn = 5;
sPreset.sPwParams.regularize_epsilon = 0.01;
sPreset.sPwParams.order = 100;
end
