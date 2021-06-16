function sDataset = GenerateDataset(actualDataDist, dim, nComponents, nTrain, nTest, interpMethod, sDatasetParams)

%% Set defaults
if ~exist('nTrain', 'var')
    nTrain = 4000;
end
if ~exist('nTest', 'var')
    nTest = 1000;
end
%% Generate data
if strcmp(actualDataDist, 'TwoMoons')
    if ~exist('sDatasetParams', 'var') || (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams, 'b_loadTwoMoonsMatFile'))
        sDatasetParams.b_loadTwoMoonsMatFile = false;
    end
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData = GenerateTwoMoonsDataset(nTrain, nTest, sDatasetParams.b_loadTwoMoonsMatFile);
    elseif strcmp(interpMethod, 'AddPoints')
        sData = GenerateTwoMoonsDataset(nTest, 0, sDatasetParams.b_loadTwoMoonsMatFile);
        data = sData.x;
        rperm = randperm(nTest);
        dataRearranged = data(rperm,:);
        sDataset.sData.x = dataRearranged(1:nTrain,:);
        sDataset.sData.xt = dataRearranged;
        labels = sData.y;
        labelsRearranged = labels(rperm,:);
        sDataset.sData.y = labelsRearranged(1:nTrain,:);
        sDataset.sData.yt = labelsRearranged;
        
    end
    omega = 0.3;
    omegaTilde = 0.3;
    dim = 2;
elseif strcmp(actualDataDist, 'TwoSpirals')
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData = GenerateTwoSpiralsDataset(nTrain, nTest);
    elseif strcmp(interpMethod, 'AddPoints')
        sData = GenerateTwoSpiralsDataset(nTest);
        data = sData.x;
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    omega = 0.3;
    omegaTilde = 0.3;
    dim = 2;
elseif strcmp(actualDataDist, 'SwissRoll')
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateSwissRoll(nTrain);
        sDataset.sData.xt = GenerateSwissRoll(nTest);
    elseif strcmp(interpMethod, 'AddPoints')
        data = GenerateSwissRoll(nTest);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    omegaTilde = 0.15;
    dim = 3;
    sDatasetParams = [];
elseif strcmp(actualDataDist, 'Gaussian')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'mu') && ~isfield(sDatasetParams,'sigma'))
        sDatasetParams.mu = 0*ones(1,dim);
        sDatasetParams.sigma = 1*eye(dim);
    end
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateGaussianData(dim, nComponents, nTrain, sDatasetParams.mu, sDatasetParams.sigma);
        sDataset.sData.xt = GenerateGaussianData(dim, nComponents, nTest, sDatasetParams.mu, sDatasetParams.sigma);
    elseif strcmp(interpMethod, 'AddPoints')
        data = GenerateGaussianData(dim, nComponents, nTest, sDatasetParams.mu, sDatasetParams.sigma);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.3*sDatasetParams.sigma(1,1);
    omegaTilde = 0.3;
elseif strcmp(actualDataDist, 'Uniform')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'xMin') && ~isfield(sDatasetParams,'xMax'))
        sDatasetParams.xMin = -1*ones(dim,1);
        sDatasetParams.xMax = 1*ones(dim,1);
    end
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateUniformData(dim, nTrain, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.xt = GenerateUniformData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
    elseif strcmp(interpMethod, 'AddPoints')
        data = GenerateUniformData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    omegaTilde = 0.15;
elseif strcmp(actualDataDist, 'Grid')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'xMin') && ~isfield(sDatasetParams,'xMax'))
        sDatasetParams.xMin = -1*ones(dim,1);
        sDatasetParams.xMax = 1*ones(dim,1);
    end
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateGridData(dim, nTrain, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.xt = GenerateGridData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
    elseif strcmp(interpMethod, 'AddPoints')
        data = GenerateGridData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
        nTest = length(data);
        r = ceil(nTest/nTrain);
        nTrain = length(1:r:nTest);
        dataRearranged(1:nTrain,:) = data(1:r:nTest,:);
        data(1:r:nTest,:) = [];
        dataRearranged(nTrain+1:nTrain+length(data),:) = data;
        
        sDataset.sData.x = dataRearranged(1:nTrain,:);
        sDataset.sData.xt = dataRearranged;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    omegaTilde = 0.15;
else
    error('unknown pdf')
end

sDataset.recommendedOmega = omega;
sDataset.recommendedOmegaTilde = omegaTilde;
sDataset.xMax = max(sDataset.sData.x);
sDataset.xMin = min(sDataset.sData.x);
sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualNumComponents = nComponents;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
sDataset.sDatasetParams = sDatasetParams;
end
