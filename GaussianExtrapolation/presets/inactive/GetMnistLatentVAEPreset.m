function sPreset = GetMnistLatentVAEPreset()
%% Dataset parameters
sPreset.dim                = 2;
sPreset.n                  = 2000;
sPreset.N                  = 2000;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'MnistLatentVAE'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'RandomWalk'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sDatasetParams.vPossibleLabels = [0, 1];
sPreset.sDatasetParams     = sDatasetParams;
%% Gaussian kernel width
sPreset.omega              = 0.1; % for nystrom kernel
sPreset.omegaTilde         = 0.1; % for our method
%% Number of signals
sPreset.nSignals           = 1;
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = length(sDatasetParams.vPossibleLabels);
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 10*sPreset.gmmNumComponents;
sPreset.MTilde             = 10*sPreset.gmmNumComponents;
%% Regularizations
sPreset.gamma1             = 0;
sPreset.gamma2             = 0;0.1;
%% Representer theorem
sPreset.gamma1Rep          = 0;
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
sPreset.b_maskDataFitTerm  = true;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'OneDraw';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end
