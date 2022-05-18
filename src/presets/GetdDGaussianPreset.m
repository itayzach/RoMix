function sPreset = GetdDGaussianPreset()
%% Dataset parameters
sPreset.dim                = 20;
sPreset.n                  = 6000;
sPreset.N                  = 6000;
sPreset.k                  = 3;
sPreset.nGenDataCompnts    = 25;
sPreset.nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
sPreset.verticesPDF        = 'Gaussian'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
sPreset.adjacencyType      = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
sPreset.matrixForEigs      = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% DatasetParams
for c = 1:sPreset.nGenDataCompnts
    sDatasetParams.mu{c}    = 5*(c-1)*ones(1,sPreset.dim);
    sDatasetParams.sigma{c} = 1*eye(sPreset.dim);
    sDatasetParams.p{c}     = 1/sPreset.nGenDataCompnts;
end
sPreset.sDatasetParams     = sDatasetParams;
%% Number of signals
sPreset.nSignals           = 1;
%% Gaussian kernel width
sPreset.omega              = sPreset.dim*0.3; % for nystrom kernel
sPreset.omegaTilde         = sPreset.dim*0.3; % for our method
%% GMM params
sPreset.gmmRegVal          = 1e-3;
sPreset.gmmMaxIter         = 2000;
sPreset.gmmNumComponents   = sPreset.nGenDataCompnts;
%% Number of eigenvectors/eigenfunctions
sPreset.M                  = 20*sPreset.gmmNumComponents;
sPreset.MTilde             = 20*sPreset.gmmNumComponents;
%% Regularizations
sPreset.gamma1             = 0.001;
sPreset.gamma2             = 0;
%% Representer theorem
sPreset.gamma1Rep          = 0.001;
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
sPreset.b_maskDataFitTerm  = false;
sPreset.b_compareMethods   = false;
%% 
sPreset.dataGenTechnique = 'OneDraw';
sPreset.sDistanceParams.distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end
