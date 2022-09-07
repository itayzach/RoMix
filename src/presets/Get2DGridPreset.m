function sPreset = Get2DGridPreset()
%% DatasetParams
sDatasetParams.xMin        = [0, 0];
sDatasetParams.xMax        = [4, 1];
sDatasetParams.nx          = [50, 25];
sDatasetParams.Nx          = [50*2, 25*2];
sPreset.sDatasetParams     = sDatasetParams;
%% Dataset parameters
sPreset.dim                = 2;
sPreset.n                  = prod(sDatasetParams.nx);
sPreset.N                  = prod(sDatasetParams.Nx);
sPreset.k                  = 3;
sPreset.nGenDataCompnts    = 0;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'Grid'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'NormLap'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% Number of signals
sPreset.nSignals           = 1;
%% Gaussian kernel width
L = norm(sDatasetParams.xMax - sDatasetParams.xMin);
sPreset.omega              = 2*sqrt(L/sPreset.n); % for nystrom kernel
sPreset.omegaTilde         = 2*sqrt(L/sPreset.n); % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = 50;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20;
sPreset.MTilde             = 20*sPreset.gmmNumComponents;1000;
%% Regularizations
sPreset.gamma1             = 1e-5;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 1e-5;
sPreset.gamma2Rep          = 0;
%% Number of runs (=realizations)
sPreset.R                  = 1;
%% Method parameters
sPreset.b_debugUseAnalytic = true;
sPreset.b_forceCtoIdentity = false;
sPreset.b_takeEigsFromWRef = false;
sPreset.b_flipSign         = true;
sPreset.b_pairwiseFlipSign = true;
sPreset.b_interpEigenvecs  = true;
sPreset.b_runGraphSignals  = true;
sPreset.b_maskDataFitTerm  = false;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'OneDraw';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end
