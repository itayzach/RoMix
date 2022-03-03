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
    
    trainTestRatio = n/N;
    N = length(data);
    n = min(round(trainTestRatio*N),N);
    N = length(data);
    rperm = randperm(N);
    dataRearranged = data(rperm,:);
    
    sDataset.sData.x = dataRearranged(1:n,:);
    sDataset.sData.xt = dataRearranged;
    
    nMonths = numel(sDatasetParams.monthNames);
    mSignals = zeros(N, nMonths);
    for monthId = 1:nMonths
        currMonth = sDatasetParams.monthNames{monthId};
        mSignals(:,monthId) = T.(currMonth)(rperm);
    end
    
    rpermUnlabeled = randperm(n);
    rpermUnlabeled = rpermUnlabeled(1:n-sDatasetParams.nLabeled);
    mSignalsMasked = mSignals(1:n,:);
    mSignalsMasked(rpermUnlabeled,:) = 0;

    sDataset.sData.y = mSignals(1:n,:);
    sDataset.sData.ymasked = mSignalsMasked;
    sDataset.sData.yt = mSignals;

elseif ismember(actualDataDist, {'USPS', 'MNIST'})
    assert(strcmp(dataGenTechnique, 'OneDraw'))
    nSets = ceil(N/1000);
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
        S_opt_prev = false(N,1);
        num_queries_to_add = sDatasetParams.nLabeled;
    
        Ln_k = 0.5*(Ln_k+Ln_k.');
        [S_opt, cutoff] = compute_opt_set_inc(Ln_k, k, num_queries_to_add, S_opt_prev);
        S_compl = setdiff((1:N)', S_opt);
        
        rperm = [S_opt.'; S_compl];
    else
%         rperm = 1:nTest;
        rperm = randperm(N);    
    end
    
    % By this point, X is ordered. We take a random permutation to randomize the data order
    Xdiv255perm = Xdiv255(rperm,:);
    mSignalsPerm = digitsDataset.mem_fn(rperm,:);
    [~, vLabels] = max(mSignalsPerm,[],2);
    
    % Take even number of samples from each class for training.
    % (There are 10 classes)
    vTrainInd = zeros(n,1);
    for classInd = 1:10
        currClassOffset = n/10*(classInd-1);
        vIndOfCurrClass = find(vLabels == classInd,n/10);
        vTrainInd(currClassOffset + (1:n/10)) = vIndOfCurrClass;
    end
    % None training indices are the complement of vTrainInt
    vNonTrainInd = setdiff((1:N)',vTrainInd);

    % Sort data as [train; test]
    Xdiv255ordered(1:n,:) = Xdiv255perm(vTrainInd,:);
    Xdiv255ordered(n+1:N,:) = Xdiv255perm(vNonTrainInd,:);
    mSignalsPermOrdered(1:n,:) = mSignalsPerm(vTrainInd,:);
    mSignalsPermOrdered(n+1:N,:) = mSignalsPerm(vNonTrainInd,:);

    % Zero out labeled signals for ssl
    nLabeled = sDatasetParams.nLabeled;
    mSignalsPermOrderedMasked = mSignalsPermOrdered(1:n,:);
    for classInd = 1:10
        vNonLabledInd = n/10*(classInd-1) + (nLabeled/10+1:n/10);
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


    sDataset.sData.x = Xdiv255ordered(1:n,:);
    sDataset.sData.xt = Xdiv255ordered;
    sDataset.sData.y = mSignalsPermOrdered(1:n,:);
    sDataset.sData.ymasked = mSignalsPermOrderedMasked;
    sDataset.sData.yt = mSignalsPermOrdered;

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
    data = double(z(1:N,:));
    labels = y(1:N)';
    sDataset.sData.x = data(1:n,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = double(labels(1:n)');
    sDataset.sData.yt = double(labels');

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
