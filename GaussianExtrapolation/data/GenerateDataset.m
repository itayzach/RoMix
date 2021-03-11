function sDataset = GenerateDataset(actualDataDist, dim, nComponents, nTrain, nTest, interpMethod, b_loadTwoMoonsMatFile)

%% Set defaults
if ~exist('b_loadTwoMoonsMatFile', 'var')
    b_loadTwoMoonsMatFile = false;
end
if ~exist('nTrain', 'var')
    nTrain = 4000;
end
if ~exist('nTest', 'var')
    nTest = 1000;
end
%% Generate data
if strcmp(actualDataDist, 'TwoMoons')
    sDataset.sData = GenerateTwoMoonsDataset(nTrain, nTest, b_loadTwoMoonsMatFile);
    omega = 0.3;
    dim = 2;
elseif strcmp(actualDataDist, 'TwoSpirals')
    sDataset.sData = GenerateTwoSpiralsDataset(nTrain, nTest);
    omega = 0.3;
    dim = 2;
elseif strcmp(actualDataDist, 'SwissRoll')
    sDataset.sData.x = GenerateSwissRoll(nTrain);
    sDataset.sData.xt = GenerateSwissRoll(nTest);
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    dim = 3;
elseif strcmp(actualDataDist, 'Gaussian')
    sDataset.sData.x = GenerateGaussianData(dim, nComponents, nTrain);
    sDataset.sData.y = [];
    sDataset.sData.xt = GenerateGaussianData(dim, nComponents, nTest);
    sDataset.sData.yt = [];
    omega = 0.3;
elseif strcmp(actualDataDist, 'Uniform')
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateUniformData(dim, nTrain);
        sDataset.sData.xt = GenerateUniformData(dim, nTest);
    elseif strcmp(interpMethod, 'AddPoints')
        data = GenerateUniformData(dim, nTest);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
elseif strcmp(actualDataDist, 'Grid')
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateGridData(dim, nTrain);
        sDataset.sData.xt = GenerateGridData(dim, nTest);
    elseif strcmp(interpMethod, 'AddPoints')
        r = round(nTest/nTrain);
        data = GenerateGridData(dim, nTest);
        dataRearranged(1:nTrain,:) = data(1:r:nTest,:);
        data(1:r:nTest,:) = [];
        dataRearranged(nTrain+1:nTest,:) = data;
        
        sDataset.sData.x = dataRearranged(1:nTrain,:);
        sDataset.sData.xt = dataRearranged;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
else
    error('unknown pdf')
end

sDataset.recommendedOmega = omega;
sDataset.xMax = max(sDataset.sData.x);
sDataset.xMin = min(sDataset.sData.x);
sDataset.nTrain = nTrain;
sDataset.nTest = nTest;
sDataset.nTotal = nTrain + nTest;
sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualNumComponents = nComponents;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
end
