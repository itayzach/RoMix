function sPreset = GetMnistCommonPreset(b_runVAE)
%% VAE parameters
latentDim                  = 20;
epochs                     = 600;
batchSize                  = 128;
nTrainVAE                  = 60000;
nTestVAE                   = 10000;
b_loadKeras                = true;
b_forceLoadTrainedVAE      = true; % Forces load of a trained model on entire MNIST (60,000 train)
%% Dataset parameters
vPer                       = [0.01 0.05 (0.1:0.2:0.7)];
sPreset.dim                = b_runVAE*latentDim + (1-b_runVAE)*28*28;
sPreset.n                  = 24000; % cannot be less than 28*28 due to GMM
sPreset.nLabeled           = round(vPer*sPreset.n);
sPreset.N                  = sPreset.n + 10000;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'MNIST'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.b_runVAE    = b_runVAE;
sDatasetParams.latentDim   = latentDim;
sDatasetParams.epochs      = epochs;
sDatasetParams.batchSize   = batchSize;
sDatasetParams.nTrainVAE   = nTrainVAE;
sDatasetParams.nTestVAE    = nTestVAE;
sDatasetParams.b_loadKeras = b_loadKeras;
sDatasetParams.b_forceLoadTrainedVAE = b_forceLoadTrainedVAE;
sDatasetParams.xTickNames  = strcat(cellfun(@num2str,(num2cell(100*vPer,numel(sPreset.nLabeled))), 'UniformOutput', false),'\%');
sDatasetParams.b_zoomInAcc = false;
sPreset.sDatasetParams     = sDatasetParams;
assert(all(sPreset.nLabeled <= sPreset.n))
%% Number of signals
sPreset.nSignals           = 1; % After conversion from number of classes
%% Gaussian kernel width
sPreset.omega              = 0.114*sPreset.dim; % for nystrom kernel
sPreset.omegaTilde         = 0.114*sPreset.dim; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-4;
sPreset.gmmMaxIter         = 200;
sPreset.gmmNumComponents   = b_runVAE*10 + (1-b_runVAE)*300;
if b_runVAE
    sPreset.vGmmNumCompAicBic  = [3, 5, 10, 15, 20];
end
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = b_runVAE*7500 + (1-b_runVAE)*1000;
%% Regularizations
sPreset.gamma1             = 0.1;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.1;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 5;
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
sPreset.cMethods = {'RoMix', 'Rep. Thm.', 'VSPW', 'Nystr\"{o}m', 'w-kNN'};
sPreset.knn = 5;
sPreset.sPwParams.regularize_epsilon = 0.01;
sPreset.sPwParams.order = 100;
end
