function sDataset = GenerateDataset(actualDataDist, dim, nTrain, nTest)

%% Set defaults
if ~exist('nTrain', 'var')
    nTrain = 4000;
end
if ~exist('nTest', 'var')
    nTest = 1000;
end
nTotal = nTrain + nTest;
%% Generate data
if strcmp(actualDataDist, 'Two_moons')
    % sDataset.sData = load('2moons.mat');
    sDataset.sData = GenerateTwoMoonsDataset(nTrain, nTest);
else
    if strcmp(actualDataDist, 'Gaussian')
        if dim == 2
            cov = [0.25    0.01; 
                   0.01   0.25];
            mu  = [0 0];
            xTotal = mvnrnd(mu, cov, nTotal);
        elseif dim == 1
            sigma = 0.5;
            mu = 100;
            
%             x = linspace(mu-5*sigma,mu+5*sigma, nTotal+1)';
%             [f, pos] = ecdf(x);
%             xTotal = icdf('Normal',f, mu, sigma);
%             xTotal = xTotal(2:end-1); % first and last samples are -Inf and Inf after icdf
            
            xTotal = sigma*randn(nTotal, 1) + mu;
        else
            error('generate random data for more than 2D')
        end        
    elseif strcmp(actualDataDist, 'Uniform')
        xTotal = (sSimParams.xMax - sSimParams.xMin)*rand(nTotal, dim) + sSimParams.xMin;
    else
        error('unknown pdf')
    end
    sDataset.sData.x = xTotal(1:nTrain,:);
    sDataset.sData.y = [];
    sDataset.sData.xt = xTotal(nTrain+1:end,:);
    sDataset.sData.yt = [];

end

sDataset.nTrain = nTrain;
sDataset.nTest = nTest;
sDataset.nTotal = nTrain + nTest;
sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
end