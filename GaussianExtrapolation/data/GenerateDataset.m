function sDataset = GenerateDataset(actualDataDist, dim, nComponents, nTrain, nTest, dataGenTechnique, sDatasetParams)
assert( (strcmp(dataGenTechnique, 'AddPoints') && nTrain <= nTest) || ...
       ~(strcmp(dataGenTechnique, 'AddPoints')))
%% Set defaults
if ~exist('nTrain', 'var')
    nTrain = 4000;
end
if ~exist('nTest', 'var')
    nTest = 1000;
end
fprintf('Generating n = %d, N = %d points of %d-d %s... ', nTrain, nTest, dim, actualDataDist)
%% Generate data
if strcmp(actualDataDist, 'SineCosine')
%     assert(strcmp(dataGenTechnique, 'AddPoints'))
%     [data, theta] = GenerateSineCosineData(dim, nTest);
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     sDataset.sData.theta = theta;
%     omega = 0.5*sum(sqrt(std(sDataset.sData.xt)))/dim; %0.15;
%     omegaTilde = 0.5*sum(sqrt(std(sDataset.sData.xt)))/dim; %0.15;
elseif strcmp(actualDataDist, 'BrazilWeather') 
    T = readtable(fullfile('data','Brazilian_Weather_Stations-Temperature_1961-1990.xlsx'));
    
    % Remove entries with NaN
    T(isnan(T.Ann),:) = [];
    
    nPoints = height(T);
    latStr = cell2mat(T.Latitude);
    lonStr = cell2mat(T.Longitude);
    % Add seconds to lat/lon to have a valid str2angle input format:
    %     '07°38''S' --> '07°38''00"S'
    latStrFormatted = [latStr(:,1:5), repmat('''00"',nPoints,1), latStr(:,7)];
    latitude = str2angle(latStrFormatted);
    lonStrFormatted = [lonStr(:,1:5), repmat('''00"',nPoints,1), lonStr(:,7)];
    longitude = str2angle(lonStrFormatted);
    data = [longitude, latitude];
    
    trainTestRatio = nTrain/nTest;
    nTest = length(data);
    nTrain = min(round(trainTestRatio*nTest),nTest);
    nTest = length(data);
    rperm = randperm(nTest);
    dataRearranged = data(rperm,:);
    
    sDataset.sData.x = dataRearranged(1:nTrain,:);
    sDataset.sData.xt = dataRearranged;
    
    nMonths = numel(sDatasetParams.monthNames);
    mSignals = zeros(nTest, nMonths);
    for monthId = 1:nMonths
        currMonth = sDatasetParams.monthNames{monthId};
        mSignals(:,monthId) = T.(currMonth)(rperm);
    end
    
    rpermUnlabeled = randperm(nTrain);
    rpermUnlabeled = rpermUnlabeled(1:nTrain-sDatasetParams.nLabeled);
    mSignalsMasked = mSignals(1:nTrain,:);
    mSignalsMasked(rpermUnlabeled,:) = 0;

    sDataset.sData.y = mSignals(1:nTrain,:);
    sDataset.sData.ymasked = mSignalsMasked;
    sDataset.sData.yt = mSignals;

elseif ismember(actualDataDist, {'USPS', 'MNIST'})
    assert(strcmp(dataGenTechnique, 'AddPoints'))
    nSets = ceil(nTest/1000);
    datasetInd = randi(10,nSets,1);
    fprintf('Loading %s set %d... ',actualDataDist, datasetInd(1))
    digitsDataset = load(['data', filesep, lower(actualDataDist), filesep, lower(actualDataDist), '_set', num2str(datasetInd(1)), '.mat']);

    for setInd = 2:nSets
        fprintf('Loading %s set %d... ',actualDataDist, datasetInd(setInd))
        digitsDataset2 = load(['data', filesep, lower(actualDataDist), filesep, lower(actualDataDist), '_set', num2str(datasetInd(setInd)), '.mat']);
        digitsDataset.X = [digitsDataset.X; digitsDataset2.X];
        digitsDataset.mem_fn = [digitsDataset.mem_fn; digitsDataset2.mem_fn];
    end
    
    Xdiv255 = digitsDataset.X/255;
    if 0
        % Set selection from "Active semisupervised learning"
        k = 8;
        sPreset = GetUspsPreset();
        [~, ~, ~, Ln] = SimpleCalcAdjacency(Xdiv255,'GaussianKernel',sPreset.sDistanceParams, sPreset.omega);
        Ln_k = Ln;
        for classInd = 1:(k-1)
            Ln_k = Ln_k*Ln;
        end
        S_opt_prev = false(nTest,1);
        num_queries_to_add = sDatasetParams.nLabeled;
    
        Ln_k = 0.5*(Ln_k+Ln_k.');
        [S_opt, cutoff] = compute_opt_set_inc(Ln_k, k, num_queries_to_add, S_opt_prev);
        S_compl = setdiff((1:nTest)', S_opt);
        
        rperm = [S_opt.'; S_compl];
    else
%         rperm = 1:nTest;
        rperm = randperm(nTest);    
    end
    
    % By this point, X is ordered. We take a random permutation to randomize the data order
    Xdiv255perm = Xdiv255(rperm,:);
    mSignalsPerm = digitsDataset.mem_fn(rperm,:);
    [~, vLabels] = max(mSignalsPerm,[],2);
    
    % Take even number of samples from each class for training.
    % (There are 10 classes)
    vTrainInd = zeros(nTrain,1);
    for classInd = 1:10
        currClassOffset = nTrain/10*(classInd-1);
        vIndOfCurrClass = find(vLabels == classInd,nTrain/10);
        vTrainInd(currClassOffset + (1:nTrain/10)) = vIndOfCurrClass;
    end
    % None training indices are the complement of vTrainInt
    vNonTrainInd = setdiff((1:nTest)',vTrainInd);

    % Sort data as [train; test]
    Xdiv255ordered(1:nTrain,:) = Xdiv255perm(vTrainInd,:);
    Xdiv255ordered(nTrain+1:nTest,:) = Xdiv255perm(vNonTrainInd,:);
    mSignalsPermOrdered(1:nTrain,:) = mSignalsPerm(vTrainInd,:);
    mSignalsPermOrdered(nTrain+1:nTest,:) = mSignalsPerm(vNonTrainInd,:);

    % Zero out labeled signals for ssl
    nLabeled = sDatasetParams.nLabeled;
    mSignalsPermOrderedMasked = mSignalsPermOrdered(1:nTrain,:);
    for classInd = 1:10
        vNonLabledInd = nTrain/10*(classInd-1) + (nLabeled/10+1:nTrain/10);
        mSignalsPermOrderedMasked(vNonLabledInd,:) = 0;
    end


%     Xdiv255_tsne = tsne(Xdiv255,'Algorithm','barneshut','NumPCAComponents',50,'NumDimensions',3);
%     [~, labels] = max(usps.mem_fn,[],2);
%     figure
%     scatter3(Xdiv255_tsne(:,1),Xdiv255_tsne(:,2),Xdiv255_tsne(:,3),15,labels,'filled')
%     figure; 
%     tiledlayout('flow')
%     imsize = sqrt(size(Xdiv255ordered,2));
%     for i = (1:30)+30*2
%         nexttile;
%         imagesc(reshape(Xdiv255ordered(i,:),imsize,imsize))
%         title(vLabelsOrdered(i))
%     end
%     sgtitle('random samples')   


    sDataset.sData.x = Xdiv255ordered(1:nTrain,:);
    sDataset.sData.xt = Xdiv255ordered;
    sDataset.sData.y = mSignalsPermOrdered(1:nTrain,:);
    sDataset.sData.ymasked = mSignalsPermOrderedMasked;
    sDataset.sData.yt = mSignalsPermOrdered;

elseif strcmp(actualDataDist, 'MnistDist')
%     assert(strcmp(dataGenTechnique, 'AddPoints'))
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
    assert(strcmp(dataGenTechnique, 'AddPoints'))
    error('change to one-hot encoding and then take argmax, like usps')
    mnistLatent = load(['data', filesep, 'mnist', filesep, 'mnistLatentVAE.mat']);
    z = double(mnistLatent.z);
    y = double(mnistLatent.labels');
    if ~isfield(sDatasetParams, 'vPossibleLabels')
        vPossibleLabels = 0:9;
    else
        vPossibleLabels = sDatasetParams.vPossibleLabels;
    end
    ind = find(ismember(y, vPossibleLabels));
    z = z(ind,:);
    y = y(ind);
    y = y + 1; % +1 since interpolating "0" is hard :)
    data = double(z(1:nTest,:));
    labels = y(1:nTest)';
    sDataset.sData.x = data(1:nTrain,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = double(labels(1:nTrain)');
    sDataset.sData.yt = double(labels');

elseif strcmp(actualDataDist, 'CoraLatentVGAE')
%     assert(strcmp(dataGenTechnique, 'AddPoints'))
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
%     assert(strcmp(dataGenTechnique, 'AddPoints'))
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
%     assert(strcmp(dataGenTechnique, 'AddPoints'))
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
%     assert(strcmp(dataGenTechnique, 'AddPoints'))
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
    if strcmp(dataGenTechnique, 'NewPoints')
        sDataset.sData = GenerateTwoMoonsDataset(nTrain, nTest, sDatasetParams.nLabeled, sDatasetParams.b_loadTwoMoonsMatFile);
    elseif strcmp(dataGenTechnique, 'AddPoints')
        sData = GenerateTwoMoonsDataset(nTest, 0, sDatasetParams.nLabeled, sDatasetParams.b_loadTwoMoonsMatFile);
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
elseif strcmp(actualDataDist, 'TwoSpirals')
    if strcmp(dataGenTechnique, 'NewPoints')
        sDataset.sData = GenerateTwoSpiralsDataset(nTrain, nTest);
    elseif strcmp(dataGenTechnique, 'AddPoints')
        sData = GenerateTwoSpiralsDataset(nTest, 0);
        data = sData.x;
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
elseif strcmp(actualDataDist, 'SwissRoll')
    if strcmp(dataGenTechnique, 'NewPoints')
        sDataset.sData.x = GenerateSwissRoll(nTrain);
        sDataset.sData.xt = GenerateSwissRoll(nTest);
    elseif strcmp(dataGenTechnique, 'AddPoints')
        data = GenerateSwissRoll(nTest);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    figure; scatter3(sDataset.sData.x(:,1),sDataset.sData.x(:,2),sDataset.sData.x(:,3),[],sDataset.sData.y,'filled');
elseif strcmp(actualDataDist, 'Gaussian')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'mu') && ~isfield(sDatasetParams,'sigma'))
        warning('mu sigma were not specified, setting defaults');
        sDatasetParams.mu = 0*ones(1,dim);
        sDatasetParams.sigma = 1*eye(dim);
    end
    if strcmp(dataGenTechnique, 'NewPoints')
        sDataset.sData.x = GenerateGaussianData(dim, nComponents, nTrain, sDatasetParams.mu, sDatasetParams.sigma);
        sDataset.sData.xt = GenerateGaussianData(dim, nComponents, nTest, sDatasetParams.mu, sDatasetParams.sigma);
    elseif strcmp(dataGenTechnique, 'AddPoints')
        data = GenerateGaussianData(dim, nComponents, nTest, sDatasetParams.mu, sDatasetParams.sigma);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end

    sDataset.sData.y = GenerateSyntheticGraphSignal(sDataset.sData.x);
    sDataset.sData.yt = GenerateSyntheticGraphSignal(sDataset.sData.xt);
    
elseif strcmp(actualDataDist, 'Uniform')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'xMin') && ~isfield(sDatasetParams,'xMax'))
        sDatasetParams.xMin = -1*ones(dim,1);
        sDatasetParams.xMax = 1*ones(dim,1);
    end
    if strcmp(dataGenTechnique, 'NewPoints')
        sDataset.sData.x = GenerateUniformData(dim, nTrain, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.xt = GenerateUniformData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
    elseif strcmp(dataGenTechnique, 'AddPoints')
        data = GenerateUniformData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.x = data(1:nTrain,:);
        sDataset.sData.xt = data;
    end
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
    
elseif strcmp(actualDataDist, 'Grid')
    if ~exist('sDatasetParams', 'var') || ...
            (exist('sDatasetParams', 'var') && ~isfield(sDatasetParams,'xMin') && ~isfield(sDatasetParams,'xMax'))
        warning('xMin xMax were not specified, setting defaults.')
        sDatasetParams.xMin = -1*ones(dim,1);
        sDatasetParams.xMax = 1*ones(dim,1);
    end
    if strcmp(dataGenTechnique, 'NewPoints')
        sDataset.sData.x = GenerateGridData(dim, nTrain, sDatasetParams.xMin, sDatasetParams.xMax);
        sDataset.sData.xt = GenerateGridData(dim, nTest, sDatasetParams.xMin, sDatasetParams.xMax);
    elseif strcmp(dataGenTechnique, 'AddPoints')
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
        
        xTrainInd = 1:r:nTest;
        sDataset.sData.x = data(xTrainInd,:);
        sDataset.sData.xt = [data(xTrainInd,:); data(setdiff(1:end,xTrainInd),:)];
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    figure; plot(sDataset.sData.x,sDataset.sData.y,'.');

    %%
    if ~isempty(msgBoxMsg)
        msgBoxTitle = 'GenerateDataset warning';
        msgBoxMsg = ['This dataset has a fixed number of points.', newline, msgBoxMsg];
        MyMsgBox(msgBoxMsg, msgBoxTitle)
    end
else
    error('unknown pdf')
end

sDataset.dataDist = actualDataDist;
sDataset.dim = dim;
sDataset.actualNumComponents = nComponents;
sDataset.actualDataDist = actualDataDist;
sDataset.estDataDist = 'Gaussian';
sDataset.sDatasetParams = sDatasetParams;

% sanity checks
assert(nTrain == size(sDataset.sData.x,1) && nTest == size(sDataset.sData.xt,1));
if strcmp(dataGenTechnique, 'NewPoints')
    warning('sPreset.dataGenTechnique is NewPoints, should make sure its okay...')
    pause(0.5)
else
    assert(isequal(sDataset.sData.x, sDataset.sData.xt(1:nTrain,:)));
end

fprintf('Done.\n')
end
