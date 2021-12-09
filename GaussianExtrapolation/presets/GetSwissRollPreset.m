function sPreset = GetSwissRollPreset()
%% Dataset parameters
sPreset.dim                = 3;
sPreset.n                  = 4000;
sPreset.N                  = 5000;
sPreset.k                  = round(0.01*sPreset.N);
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'SwissRoll'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'RandomWalk'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
sPreset.sDatasetParams     = [];
%% Gaussian kernel width
sPreset.omega              = sqrt(5)/sqrt(2); % for nystrom kernel
sPreset.omegaTilde         = sqrt(5)/sqrt(2); % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 20;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 9;
sPreset.MTilde             = 25*sPreset.gmmNumComponents;
sPreset.gamma1             = 0; %5;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0;
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