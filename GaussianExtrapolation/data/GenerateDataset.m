function sDataset = GenerateDataset(actualDataDist, dim, nComponents, nTrain, nTest, b_loadTwoMoonsMatFile)

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
nTotal = nTrain + nTest;
%% Generate data
if strcmp(actualDataDist, 'TwoMoons')
    sDataset.sData = GenerateTwoMoonsDataset(nTrain, nTest, b_loadTwoMoonsMatFile);
    omega = 0.3;
    dim = 2;
    xMax = max(sDataset.sData.x);
    xMin = min(sDataset.sData.x);
elseif strcmp(actualDataDist, 'TwoSpirals')
    sDataset.sData = GenerateTwoSpiralsDataset(nTrain, nTest);
    omega = 0.3;
    dim = 2;
    xMax = max(sDataset.sData.x);
    xMin = min(sDataset.sData.x);
elseif strcmp(actualDataDist, 'SwissRoll')
    sDataset.sData.x = GenerateSwissRoll(nTrain);
    sDataset.sData.xt = GenerateSwissRoll(nTest);
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.3;
    dim = 3;
    xMax = max(sDataset.sData.x);
    xMin = min(sDataset.sData.x);
else
    if strcmp(actualDataDist, 'Gaussian')
        if dim == 2
            if nComponents == 1
                cov = [0.25    0.01; 
                       0.01   0.25];
                mu  = [0; 0];
                xTotal = mvnrnd(mu, cov, nTotal);
            elseif nComponents == 2
                cov1 = [0.25    0.01; 
                        0.01   0.5];
                mu1  = [0; 0];
                x1 = mvnrnd(mu1, cov1, nTotal);
                
                cov2 = [0.7    -0.1; 
                        -0.1   0.1];
                mu2  = [5; 5];
                x2 = mvnrnd(mu2, cov2, nTotal);
                
                vSel = rand(nTotal, 1) < 0.5;
                xTotal = (1-vSel).*x1 + vSel.*x2;
            else
                error('not supported')
            end
        elseif dim == 1
            if nComponents == 1
                sigma = 1;
                mu = 0;
                xTotal = sigma*randn(nTotal, 1) + mu;
                xMax = 3*sigma;
                xMin = -3*sigma;
                omega = 0.3;
            elseif nComponents == 2
                mu = [2; 7];       % Means
                sigma = [0.4 0.5]; % Covariances
                vSel = rand(nTotal, 1) < 0.5;
                xTotal = (1-vSel).*(sigma(1)*randn(nTotal, 1) + mu(1)) + ...
                           vSel.*(sigma(2)*randn(nTotal, 1) + mu(2));
            else
                error('not supported')
            end
        else
            error('generate random data for more than 2D')
        end        
    elseif strcmp(actualDataDist, 'Uniform')
        xMin = -1;
        xMax = 1;
        xTotal = (xMax - xMin)*rand(nTotal, dim) + xMin;
        omega = 0.3;
    else
        error('unknown pdf')
    end
    sDataset.sData.x = xTotal(1:nTrain,:);
    sDataset.sData.y = [];
    sDataset.sData.xt = xTotal(nTrain+1:end,:);
    sDataset.sData.yt = [];

end

sDataset.recommendedOmega = omega;
sDataset.xMax = xMax;
sDataset.xMin = xMin;
sDataset.nTrain = nTrain;
sDataset.nTest = nTest;
sDataset.nTotal = nTrain + nTest;
sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualNumComponents = nComponents;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
end
