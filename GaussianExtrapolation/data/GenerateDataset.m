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
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData = GenerateTwoMoonsDataset(nTrain, nTest, b_loadTwoMoonsMatFile);
    elseif strcmp(interpMethod, 'AddPoints')
        sData = GenerateTwoMoonsDataset(nTest, 0, b_loadTwoMoonsMatFile);
        data = sData.x;
        dataRearranged = data(randperm(nTest),:);
        sDataset.sData.x = dataRearranged(1:nTrain,:);
        sDataset.sData.xt = dataRearranged;
    end
    omega = 0.3;
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
    dim = 3;
elseif strcmp(actualDataDist, 'Gaussian')
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData.x = GenerateGaussianData(dim, nComponents, nTrain);
        sDataset.sData.xt = GenerateGaussianData(dim, nComponents, nTest);
    elseif strcmp(interpMethod, 'AddPoints')
        data = GenerateGaussianData(dim, nComponents, nTest);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    sDataset.sData.y = [];
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
else
    error('unknown pdf')
end

sDataset.recommendedOmega = omega;
sDataset.xMax = max(sDataset.sData.x);
sDataset.xMin = min(sDataset.sData.x);
sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualNumComponents = nComponents;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
end
