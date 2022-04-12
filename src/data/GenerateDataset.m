function sDataset = GenerateDataset(actualDataDist, dim, nComponents, n, N, dataGenTechnique, sDatasetParams)
assert( (strcmp(dataGenTechnique, 'OneDraw') && n <= N) || ...
       ~(strcmp(dataGenTechnique, 'OneDraw')))
%% Generate data
fprintf('Generating n = %d, N = %d points of %d-d %s... ', n, N, dim, actualDataDist)
if strcmp(actualDataDist, 'SineCosine')
%     assert(strcmp(dataGenTechnique, 'OneDraw'))
%     [data, theta] = GenerateSineCosineData(dim, nTest);
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     sDataset.sData.theta = theta;
%     omega = 0.5*sum(sqrt(std(sDataset.sData.xt)))/dim; %0.15;
%     omegaTilde = 0.5*sum(sqrt(std(sDataset.sData.xt)))/dim; %0.15;
elseif strcmp(actualDataDist, 'BrazilWeather') 
    sDataset = LoadBrazilWeatherDataset(sDatasetParams, N, n);

elseif ismember(actualDataDist, {'USPS', 'MNIST'})
    assert(strcmp(dataGenTechnique, 'OneDraw'))
    sDataset = LoadDigitsDataset(actualDataDist, sDatasetParams, N, n);
    
elseif strcmp(actualDataDist, 'MnistDist')
%     assert(strcmp(dataGenTechnique, 'OneDraw'))
%     mnist = load('data\mnist.mat');
%     dist = pdist2(mnist.testX(1:nTest,:), mnist.testX(1:nTest,:));
%     data = triu(dist,1);
%     data = data(data > 0);
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     dim = 2;
elseif strcmp(actualDataDist, 'MnistLatentVAE')
    assert(strcmp(dataGenTechnique, 'OneDraw'))
    %error('change to one-hot encoding and then take argmax, like usps')
    mnistLatent = load(['data', filesep, 'mnist', filesep, 'mnistLatentVAE', num2str(dim), 'd.mat']);
    z = double(mnistLatent.z);
    labels = double(mnistLatent.labels');
    if ~isfield(sDatasetParams, 'vPossibleLabels')
        vPossibleLabels = 0:9;
    else
        vPossibleLabels = sDatasetParams.vPossibleLabels;
    end
    ind = find(ismember(labels, vPossibleLabels));
    z = z(ind,:);
    labels = labels(ind);
    data = double(z(1:N,:));
    labels = labels(1:N);

    oneHotEnc = false(N,10);
    for i = vPossibleLabels
        oneHotEnc(labels==i,i+1) = true;
    end
    sDataset.sData.x = data(1:n,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = double(oneHotEnc(1:n,:));
    sDataset.sData.yt = double(oneHotEnc);

elseif strcmp(actualDataDist, 'CoraLatentVGAE')
%     assert(strcmp(dataGenTechnique, 'OneDraw'))
%     coraLatent = load('data\coraLatentVGAE.mat');
%     z = double(coraLatent.z_mean);
%     data = double(z(1:nTest,:));
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     omega = 1.5;
%     omegaTilde = 1.5;
%     dim = size(data, 2);
elseif strcmp(actualDataDist, 'Minnesota')   
%     assert(strcmp(dataGenTechnique, 'OneDraw'))
%     G = gsp_minnesota();
%     data = G.coords;
%     trainTestRatio = nTrain/nTest;
%     nTest = length(data);
%     nTrain = min(round(trainTestRatio*nTest),nTest);
%     
%     rperm = randperm(nTest);
%     dataRearranged = data(rperm,:);
%    
%     sDataset.sData.x = dataRearranged(1:nTrain,:);
%     sDataset.sData.xt = dataRearranged;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     omega = 0.5;
%     omegaTilde = 0.5;
%     dim = 2;
%     sDatasetParams = [];
%     msgBoxTitle = 'GenerateDataset warning';
%     msgBoxMsg = ['This dataset has a fixed number of points.', newline, ...
%         'updated to N = ' num2str(nTest), ', n = ', num2str(nTrain)];
%     MyMsgBox(msgBoxMsg, msgBoxTitle)
elseif strcmp(actualDataDist, 'Bunny')
%     assert(strcmp(dataGenTechnique, 'OneDraw'))
%     G = gsp_bunny();
%     data = G.coords;
%     trainTestRatio = nTrain/nTest;
%     nTest = length(data);
%     nTrain = min(round(trainTestRatio*nTest),nTest);
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     omega = 0.15;
%     omegaTilde = 0.15;
%     dim = 3;
%     sDatasetParams = [];
%     msgBoxTitle = 'GenerateDataset warning';
%     msgBoxMsg = ['This dataset has a fixed number of points.', newline, ...
%         'updated to N = ' num2str(nTest), ', n = ', num2str(nTrain)];
%     MyMsgBox(msgBoxMsg, msgBoxTitle)
%     msgBoxTitle = 'GenerateDataset warning';
%     msgBoxMsg = 'This dataset requires omega adjustments...';
%     MyMsgBox(msgBoxMsg, msgBoxTitle)
elseif strcmp(actualDataDist, 'Sensor')
%     assert(strcmp(dataGenTechnique, 'OneDraw'))
%     G = gsp_sensor(nTest);
%     data = G.coords;
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     omega = 0.15;
%     omegaTilde = 0.15;
%     dim = 2;
elseif strcmp(actualDataDist, 'TwoMoons')
    if ~exist('sDatasetParams', 'var') || (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams, 'b_loadTwoMoonsMatFile'))
        sDatasetParams.b_loadTwoMoonsMatFile = false;
    end
    if strcmp(dataGenTechnique, 'TwoDraws')
        sDataset.sData = GenerateTwoMoonsDataset(n, N, sDatasetParams.nLabeled, sDatasetParams.b_loadTwoMoonsMatFile);
    elseif strcmp(dataGenTechnique, 'OneDraw')
        sData = GenerateTwoMoonsDataset(N, 0, sDatasetParams.nLabeled, sDatasetParams.b_loadTwoMoonsMatFile);
        data = sData.x;
        rperm = randperm(N);
        dataRearranged = data(rperm,:);
        sDataset.sData.x = dataRearranged(1:n,:);
        sDataset.sData.xt = dataRearranged;
        labels = sData.y;
        labelsRearranged = labels(rperm,:);
        sDataset.sData.y = labelsRearranged(1:n,:);
        sDataset.sData.yt = labelsRearranged;
    end
elseif strcmp(actualDataDist, 'TwoSpirals')
    if strcmp(dataGenTechnique, 'TwoDraws')
        sDataset.sData = GenerateTwoSpiralsDataset(n, N);
    elseif strcmp(dataGenTechnique, 'OneDraw')
        sData = GenerateTwoSpiralsDataset(N, 0);
        data = sData.x;
        sDataset.sData.x = data(1:n,:);
        sDataset.sData.xt = data;
    end
elseif strcmp(actualDataDist, 'SwissRoll')
    if strcmp(dataGenTechnique, 'TwoDraws')
        [sDataset.sData.x, sDataset.sData.S] = GenerateSwissRoll(n, sDatasetParams.a, sDatasetParams.maxTheta, sDatasetParams.height);
        [sDataset.sData.xt, sDataset.sData.St] = GenerateSwissRoll(N, sDatasetParams.a, sDatasetParams.maxTheta, sDatasetParams.height);
    elseif strcmp(dataGenTechnique, 'OneDraw')
        [data, S] = GenerateSwissRoll(N, sDatasetParams.a, sDatasetParams.maxTheta, sDatasetParams.height);
        sDataset.sData.x = data(1:n,:);
        sDataset.sData.xt = data;
        sDataset.sData.S = S(1:n,:);
        sDataset.sData.St = S;
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.S, sDataset.sData.St);
    %figure; scatter3(sDataset.sData.x(:,1),sDataset.sData.x(:,2),sDataset.sData.x(:,3),[],sDataset.sData.y,'filled');
elseif strcmp(actualDataDist, 'Gaussian')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'mu') && ~isfield(sDatasetParams,'sigma'))
        warning('mu sigma were not specified, setting defaults');
        sDatasetParams.mu = 0*ones(1,dim);
        sDatasetParams.sigma = 1*eye(dim);
    end
    if strcmp(dataGenTechnique, 'TwoDraws')
        sDataset.sData.x = GenerateGaussianData(dim, nComponents, n, sDatasetParams.mu, sDatasetParams.sigma);
        sDataset.sData.xt = GenerateGaussianData(dim, nComponents, N, sDatasetParams.mu, sDatasetParams.sigma);
    elseif strcmp(dataGenTechnique, 'OneDraw')
        data = GenerateGaussianData(dim, nComponents, N, sDatasetParams.mu, sDatasetParams.sigma);
        sDataset.sData.x = data(1:n,:);
        sDataset.sData.xt = data;
    end

    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    %figure; plot(sDataset.sData.x,sDataset.sData.y,'o',sDataset.sData.xt,sDataset.sData.yt,'.');
elseif strcmp(actualDataDist, 'Uniform')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'xMin') && ~isfield(sDatasetParams,'xMax'))
        sDatasetParams.xMin = -1*ones(dim,1);
        sDatasetParams.xMax = 1*ones(dim,1);
    end
    if strcmp(dataGenTechnique, 'TwoDraws')
        sDataset.sData.x = GenerateUniformData(dim, n, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.xt = GenerateUniformData(dim, N, sDatasetParams.xMin, sDatasetParams.xMax);
    elseif strcmp(dataGenTechnique, 'OneDraw')
        data = GenerateUniformData(dim, N, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.x = data(1:n,:);
        sDataset.sData.xt = data;
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    
elseif strcmp(actualDataDist, 'Grid')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'xMin') && ~isfield(sDatasetParams,'xMax'))
        warning('xMin xMax were not specified, setting defaults.')
        sDatasetParams.xMin = -1*ones(dim,1);
        sDatasetParams.xMax = 1*ones(dim,1);
    end
    if strcmp(dataGenTechnique, 'TwoDraws')
        sDataset.sData.x = GenerateGridData(dim, sDatasetParams.nx, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.xt = GenerateGridData(dim, sDatasetParams.Nx, sDatasetParams.xMin, sDatasetParams.xMax);
    elseif strcmp(dataGenTechnique, 'OneDraw')
        data = GenerateGridData(dim, sDatasetParams.Nx, sDatasetParams.xMin, sDatasetParams.xMax);
        msgBoxMsg = [];
        if N ~= length(data)
            N = length(data);
            msgBoxMsg = [msgBoxMsg, 'updated to N = ' num2str(N)];
        end
        r = ceil(N/n);
        if n ~= length(1:r:N)
            n = length(1:r:N);
            msgBoxMsg = [msgBoxMsg, 'updated to n = ' num2str(n)];
        end
        
        if dim == 1
            xTrainInd = 1:r:N;
        elseif dim == 2
            assert(sqrt(r) == floor(sqrt(r)));
            [col, row] = meshgrid(1:sqrt(r):sDatasetParams.Nx(1),1:sqrt(r):sDatasetParams.Nx(2));
            xTrainInd = sub2ind(fliplr(sDatasetParams.Nx),row(:),col(:));
        else
            error('invalid dim for this dataset')
        end
        sDataset.sData.x = data(xTrainInd,:);
        sDataset.sData.xt = [data(xTrainInd,:); data(setdiff(1:end,xTrainInd),:)];
        if ~isempty(msgBoxMsg)
            msgBoxTitle = 'GenerateDataset warning';
            msgBoxMsg = ['This dataset has a fixed number of points.', newline, msgBoxMsg];
            MyMsgBox(msgBoxMsg, msgBoxTitle)
        end
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    %figure; plot(sDataset.sData.x(:,1),sDataset.sData.y,'.');
    %figure; scatter3(sDataset.sData.x(:,1),sDataset.sData.x(:,2),sDataset.sData.y,[],sDataset.sData.y,'filled');
    %figure; scatter(sDataset.sData.xt(:,1),sDataset.sData.xt(:,2)); hold on; scatter(sDataset.sData.x(:,1),sDataset.sData.x(:,2),'filled')

else
    error('unknown pdf')
end

sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualNumComponents = nComponents;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
sDataset.sDataParams = sDatasetParams;

% sanity checks
assert(n == size(sDataset.sData.x,1) && N == size(sDataset.sData.xt,1));
if strcmp(dataGenTechnique, 'TwoDraws')
    warning('sPreset.dataGenTechnique is TwoDraws, should make sure its okay...')
    pause(0.5)
else
    assert(isequal(sDataset.sData.x, sDataset.sData.xt(1:n,:)));
end

fprintf('Done.\n')
end
