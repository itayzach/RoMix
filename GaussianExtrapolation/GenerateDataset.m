function sDataset = GenerateDataset(actualDataDist, nTrain, nTest)

%% Set defaults
if ~exist('nTrain', 'var')
    nTrain = 4000;
end
if ~exist('nTest', 'var')
    nTest = 1000;
end
nTotal = nTrain + nTest;
dim = 2;
%% Generate data
if strcmp(actualDataDist, 'Two_moons')
    % sDataset.mData = load('2moons.mat');
    sDataset.mData = GenerateTwoMoonsDataset(nTrain, nTest);
else
    if strcmp(actualDataDist, 'Gaussian')
        if dim == 2
            cov = [0.25    0.01; 
                   0.01   0.25];
            mu  = [0 0];
        elseif dim == 1
            cov = 0.5;
            mu = 100;
        else
            error('generate random data for more than 2D')
        end        
        xTotal = mvnrnd(mu, cov, nTotal);
    elseif strcmp(actualDataDist, 'Uniform')
        xTotal = (sSimParams.xMax - sSimParams.xMin)*rand(nTotal, dim) + sSimParams.xMin;
    else
        error('unknown pdf')
    end
    sDataset.mData.x = xTotal(1:nTrain,:);
    sDataset.mData.y = [];
    sDataset.mData.xt = xTotal(nTrain+1:end,:);
    sDataset.mData.yt = [];

end

sDataset.nTrain = nTrain;
sDataset.nTest = nTest;
sDataset.nTotal = nTrain + nTest;
sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
end