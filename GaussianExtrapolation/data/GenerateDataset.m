function sDataset = GenerateDataset(actualDataDist, dim, nComponents, nTrain, nTest, interpMethod, sDatasetParams)

%% Set defaults
if ~exist('nTrain', 'var')
    nTrain = 4000;
end
if ~exist('nTest', 'var')
    nTest = 1000;
end
%% Generate data
if strcmp(actualDataDist, 'SineCosine')
    assert(strcmp(interpMethod, 'AddPoints'))
    [data, theta] = GenerateSineCosineData(dim, nTest);
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    sDataset.sData.theta = theta;
    omega = 0.5*sum(sqrt(std(sDataset.sData.xt)))/dim; %0.15;
    omegaTilde = 0.5*sum(sqrt(std(sDataset.sData.xt)))/dim; %0.15;
    
elseif strcmp(actualDataDist, 'MnistDist')
    assert(strcmp(interpMethod, 'AddPoints'))
    mnist = load('data\mnist.mat');
    dist = pdist2(mnist.testX(1:nTest,:), mnist.testX(1:nTest,:));
    data = triu(dist,1);
    data = data(data > 0);
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    dim = 2;
    sDatasetParams = [];
elseif strcmp(actualDataDist, 'MnistLatentVAE')
    assert(strcmp(interpMethod, 'AddPoints'))
    mnistLatent = load('data\mnistLatentVAE.mat');
    z = double(mnistLatent.z);
    y = double(mnistLatent.labels');
    possibleLabels = 0:9; %[0, 1, 4];
    ind = find(ismember(y, possibleLabels));
    z = z(ind,:);
    y = y(ind);
    data = double(z(1:nTest,:));
    labels = y(1:nTest)';
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = double(labels(1:nTrain)');
    sDataset.sData.yt = double(labels');
    omega = 0.5;
    omegaTilde = 0.5;
    dim = size(data, 2);
    sDatasetParams = [];
elseif strcmp(actualDataDist, 'Minnesota')   
    assert(strcmp(interpMethod, 'AddPoints'))
    G = gsp_minnesota();
    data = G.coords;
    trainTestRatio = nTrain/nTest;
    nTest = length(data);
    nTrain = min(round(trainTestRatio*nTest),nTest);
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    omegaTilde = 0.15;
    dim = 2;
    sDatasetParams = [];
    msgBoxTitle = 'GenerateDataset warning';
    msgBoxMsg = ['This dataset has a fixed number of points.', newline, ...
        'updated to N = ' num2str(nTest), ', n = ', num2str(nTrain)];
    MyMsgBox(msgBoxMsg, msgBoxTitle)
elseif strcmp(actualDataDist, 'Bunny')
    assert(strcmp(interpMethod, 'AddPoints'))
    G = gsp_bunny();
    data = G.coords;
    trainTestRatio = nTrain/nTest;
    nTest = length(data);
    nTrain = min(round(trainTestRatio*nTest),nTest);
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    omegaTilde = 0.15;
    dim = 3;
    sDatasetParams = [];
    msgBoxTitle = 'GenerateDataset warning';
    msgBoxMsg = ['This dataset has a fixed number of points.', newline, ...
        'updated to N = ' num2str(nTest), ', n = ', num2str(nTrain)];
    MyMsgBox(msgBoxMsg, msgBoxTitle)
    msgBoxTitle = 'GenerateDataset warning';
    msgBoxMsg = 'This dataset requires omega adjustments...';
    MyMsgBox(msgBoxMsg, msgBoxTitle)
elseif strcmp(actualDataDist, 'Sensor')
    assert(strcmp(interpMethod, 'AddPoints'))
    G = gsp_sensor(nTest);
    data = G.coords;
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    omega = 0.15;
    omegaTilde = 0.15;
    dim = 2;
    sDatasetParams = [];    
elseif strcmp(actualDataDist, 'TwoMoons')
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
    dim = 2;
    omega = 0.15;
    omegaTilde = 0.15;
elseif strcmp(actualDataDist, 'TwoSpirals')
    if strcmp(interpMethod, 'NewPoints')
        sDataset.sData = GenerateTwoSpiralsDataset(nTrain, nTest);
    elseif strcmp(interpMethod, 'AddPoints')
        sData = GenerateTwoSpiralsDataset(nTest, 0);
        data = sData.x;
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    dim = 2;
    omega = 0.15;
    omegaTilde = 0.15;
    msgBoxTitle = 'GenerateDataset warning';
    msgBoxMsg = 'This dataset requires omega adjustments...';
    MyMsgBox(msgBoxMsg, msgBoxTitle)
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
    dim = 3;
    omega = sqrt(5)/sqrt(2); %0.15;
    omegaTilde = sqrt(5)/sqrt(2); %0.15;
%     omega = 2*sqrt(sum(std(sDataset.sData.xt))/dim);
%     omegaTilde = 2*sqrt(sum(std(sDataset.sData.xt))/dim);
    
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
%     warning('omega = 0.3*sigma(1,1)');
%     omega = 0.3*sDatasetParams.sigma(1,1);
    omega = 0.3;
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
        msgBoxMsg = [];
        if nTest ~= length(data)
            nTest = length(data);
            msgBoxMsg = [msgBoxMsg, 'updated to N = ' num2str(nTest)];
        end
        r = ceil(nTest/nTrain);
        if nTrain ~= length(1:r:nTest)
            nTrain = length(1:r:nTest);
            msgBoxMsg = [msgBoxMsg, 'updated to n = ' num2str(nTrain)];
        end
        
        dataRearranged(1:nTrain,:) = data(1:r:nTest,:);
        data(1:r:nTest,:) = [];
        dataRearranged(nTrain+1:nTrain+length(data),:) = data;        
        sDataset.sData.x = dataRearranged(1:nTrain,:);
        sDataset.sData.xt = dataRearranged;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    dist = norm(data(2,:)-data(1,:));
    omega = 2*sqrt(dist);
    omegaTilde = 2*sqrt(dist);
%     omega = 0.5*std(sDataset.sData.xt)^2;
%     omegaTilde = 0.5*std(sDataset.sData.xt)^2;
%     omega = 0.15;
%     omegaTilde = 0.15;
    if ~isempty(msgBoxMsg)
        msgBoxTitle = 'GenerateDataset warning';
        msgBoxMsg = ['This dataset has a fixed number of points.', newline, msgBoxMsg];
        MyMsgBox(msgBoxMsg, msgBoxTitle)
    end
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
