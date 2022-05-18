function sPreset = GetMnistCommonPreset(b_runVAE)
%% VAE parameters
latentDim                  = 20;
epochs                     = 200;
batchSize                  = 128;
b_loadKeras                = true;
b_forceLoadTrainedVAE      = true; % Forces load of a trained model on entire MNIST (60,000 train)
%% Dataset parameters
sPreset.dim                = b_runVAE*latentDim + (1-b_runVAE)*28*28;
sPreset.n                  = 8000; % cannot be less than 28*28 due to GMM
sPreset.N                  = 10000;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'MNIST'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.nLabeled    = round(0.5*sPreset.n);
sDatasetParams.b_runVAE    = b_runVAE;
sDatasetParams.latentDim   = latentDim;
sDatasetParams.epochs      = epochs;
sDatasetParams.batchSize   = batchSize;
sDatasetParams.b_loadKeras = b_loadKeras;
sDatasetParams.b_forceLoadTrainedVAE = b_forceLoadTrainedVAE;
sPreset.sDatasetParams     = sDatasetParams;
assert(sDatasetParams.nLabeled <= sPreset.n)
%% Number of signals
sPreset.nSignals           = 1; % After conversion from number of classes
%% Gaussian kernel width
sPreset.omega              = 0.2*sPreset.dim; % for nystrom kernel
sPreset.omegaTilde         = 0.2*sPreset.dim; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = b_runVAE*1 + (1-b_runVAE)*300;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = b_runVAE*3000 + (1-b_runVAE)*1000;
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
end
