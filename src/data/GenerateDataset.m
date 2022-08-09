function sDataset = GenerateDataset(sPlotParams, sPreset, b_interpEigenvecs)
assert( (strcmp(sPreset.dataGenTechnique, 'OneDraw') && sPreset.n <= sPreset.N) || ...
       ~(strcmp(sPreset.dataGenTechnique, 'OneDraw')))
%% Generate data
fprintf('Generating n = %d, N = %d points of %d-d %s... ', sPreset.n, sPreset.N, sPreset.dim, sPreset.verticesPDF)
if strcmp(sPreset.verticesPDF, 'SineCosine')
%     assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
%     [data, theta] = GenerateSineCosineData(sPreset.dim, nTest);
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     sDataset.sData.theta = theta;
%     omega = 0.5*sum(sqrt(std(sDataset.sData.xt)))/sPreset.dim; %0.15;
%     omegaTilde = 0.5*sum(sqrt(std(sDataset.sData.xt)))/sPreset.dim; %0.15;
elseif strcmp(sPreset.verticesPDF, 'BrazilWeather') 
    sDataset = LoadBrazilWeatherDataset(sPreset.sDatasetParams, sPreset.N, sPreset.n, sPreset.nLabeled);

elseif ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})
    assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
    if isfield(sPreset.sDatasetParams,'b_loadKeras') && sPreset.sDatasetParams.b_loadKeras
        assert(strcmp(sPreset.verticesPDF, 'MNIST'))
        sDataset = LoadDigitsKeras(sPlotParams, sPreset.sDatasetParams, sPreset.N, sPreset.n, sPreset.nLabeled);
    else
        sDataset = LoadDigitsDataset(sPlotParams,sPreset.verticesPDF, sPreset.sDatasetParams, sPreset.N, sPreset.n);
    end
    
elseif strcmp(sPreset.verticesPDF, 'MnistDist')
%     assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
%     mnist = load('data\mnist.mat');
%     dist = pdist2(mnist.testX(1:nTest,:), mnist.testX(1:nTest,:));
%     data = triu(dist,1);
%     data = data(data > 0);
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     sPreset.dim = 2;
elseif strcmp(sPreset.verticesPDF, 'MnistLatentVAE')
    assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
    %error('change to one-hot encoding and then take argmax, like usps')
    sMnistLatent = load(['data', filesep, 'mnist', filesep, 'mnistLatentVAE', num2str(sPreset.dim), 'd.mat']);
    z = double(sMnistLatent.z);
    labels = double(sMnistLatent.labels');
    if ~isfield(sPreset.sDatasetParams, 'vPossibleLabels')
        vPossibleLabels = 0:9;
    else
        vPossibleLabels = sPreset.sDatasetParams.vPossibleLabels;
    end
    ind = find(ismember(labels, vPossibleLabels));
    z = z(ind,:);
    labels = labels(ind);
    data = double(z(1:sPreset.N,:));
    labels = labels(1:sPreset.N);

    oneHotEnc = false(sPreset.N,10);
    for i = vPossibleLabels
        oneHotEnc(labels==i,i+1) = true;
    end
    sDataset.sData.x = data(1:sPreset.n,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = double(oneHotEnc(1:sPreset.n,:));
    sDataset.sData.yt = double(oneHotEnc);

elseif strcmp(sPreset.verticesPDF, 'CoraLatentVGAE')
    assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
    sCoraLatent = load('data\coraLatentVGAE.mat');
    z = double(sCoraLatent.z_mean);
    data = double(z(1:sPreset.N,:));
    sDataset.sData.x = data(1:sPreset.n,:);
    sDataset.sData.xt = data;
    sDataset.sData.y = [];
    sDataset.sData.yt = [];
elseif strcmp(sPreset.verticesPDF, 'Minnesota')   
%     assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
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
%     sPreset.dim = 2;
%     sPreset.sDatasetParams = [];
%     msgBoxTitle = 'GenerateDataset warning';
%     msgBoxMsg = ['This dataset has a fixed number of points.', newline, ...
%         'updated to sPreset.N = ' num2str(nTest), ', sPreset.n = ', num2str(nTrain)];
%     MyMsgBox(msgBoxMsg, msgBoxTitle)
elseif strcmp(sPreset.verticesPDF, 'Bunny')
%     assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
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
%     sPreset.dim = 3;
%     sPreset.sDatasetParams = [];
%     msgBoxTitle = 'GenerateDataset warning';
%     msgBoxMsg = ['This dataset has a fixed number of points.', newline, ...
%         'updated to sPreset.N = ' num2str(nTest), ', sPreset.n = ', num2str(nTrain)];
%     MyMsgBox(msgBoxMsg, msgBoxTitle)
%     msgBoxTitle = 'GenerateDataset warning';
%     msgBoxMsg = 'This dataset requires omega adjustments...';
%     MyMsgBox(msgBoxMsg, msgBoxTitle)
elseif strcmp(sPreset.verticesPDF, 'Sensor')
%     assert(strcmp(sPreset.dataGenTechnique, 'OneDraw'))
%     G = gsp_sensor(nTest);
%     data = G.coords;
%     sDataset.sData.x = data(1:nTrain,:);
%     sDataset.sData.xt = data;
%     sDataset.sData.y = [];
%     sDataset.sData.yt = [];
%     omega = 0.15;
%     omegaTilde = 0.15;
%     sPreset.dim = 2;
elseif strcmp(sPreset.verticesPDF, 'TwoMoons')
    if ~isfield(sPreset, 'sDatasetParams') || (isfield(sPreset, 'sDatasetParams') && ~isfield(sPreset.sDatasetParams, 'b_loadTwoMoonsMatFile'))
        sPreset.sDatasetParams.b_loadTwoMoonsMatFile = false;
    end
    if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
        sDataset.sData = GenerateTwoMoonsDataset(sPreset.n, sPreset.N, sPreset.nLabeled, sPreset.sDatasetParams.b_loadTwoMoonsMatFile);
    elseif strcmp(sPreset.dataGenTechnique, 'OneDraw')
        sData = GenerateTwoMoonsDataset(sPreset.N, 0, sPreset.nLabeled, sPreset.sDatasetParams.b_loadTwoMoonsMatFile);
        data = sData.x;
        rperm = randperm(sPreset.N);
        dataRearranged = data(rperm,:);
        sDataset.sData.x = dataRearranged(1:sPreset.n,:);
        sDataset.sData.xt = dataRearranged;
        labels = sData.y;
        labelsRearranged = labels(rperm,:);
        sDataset.sData.y = labelsRearranged(1:sPreset.n,:);
        sDataset.sData.yt = labelsRearranged;
    end
elseif strcmp(sPreset.verticesPDF, 'TwoSpirals')
    if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
        sDataset.sData = GenerateTwoSpiralsDataset(sPreset.n, sPreset.N);
    elseif strcmp(sPreset.dataGenTechnique, 'OneDraw')
        sData = GenerateTwoSpiralsDataset(sPreset.N, 0);
        data = sData.x;
        sDataset.sData.x = data(1:sPreset.n,:);
        sDataset.sData.xt = data;
    end
elseif strcmp(sPreset.verticesPDF, 'SwissRoll')
    if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
        [sDataset.sData.x, sDataset.sData.S] = GenerateSwissRoll(sPreset.n, sPreset.sDatasetParams.a, sPreset.sDatasetParams.maxTheta, sPreset.sDatasetParams.height);
        [sDataset.sData.xt, sDataset.sData.St] = GenerateSwissRoll(sPreset.N, sPreset.sDatasetParams.a, sPreset.sDatasetParams.maxTheta, sPreset.sDatasetParams.height);
    elseif strcmp(sPreset.dataGenTechnique, 'OneDraw')
        [data, S] = GenerateSwissRoll(sPreset.N, sPreset.sDatasetParams.a, sPreset.sDatasetParams.maxTheta, sPreset.sDatasetParams.height);
        sDataset.sData.x = data(1:sPreset.n,:);
        sDataset.sData.xt = data;
        sDataset.sData.S = S(1:sPreset.n,:);
        sDataset.sData.St = S;
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.S, sDataset.sData.St);
    %figure; scatter3(sDataset.sData.x(:,1),sDataset.sData.x(:,2),sDataset.sData.x(:,3),[],sDataset.sData.y,'filled');
elseif strcmp(sPreset.verticesPDF, 'Gaussian')
    if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
        sDataset.sData.x = GenerateGaussianData(sPreset.dim, sPreset.nGenDataCompnts, sPreset.n, sPreset.sDatasetParams.mu, sPreset.sDatasetParams.sigma, sPreset.sDatasetParams.p);
        sDataset.sData.xt = GenerateGaussianData(sPreset.dim, sPreset.nGenDataCompnts, sPreset.N, sPreset.sDatasetParams.mu, sPreset.sDatasetParams.sigma, sPreset.sDatasetParams.p);
    elseif strcmp(sPreset.dataGenTechnique, 'OneDraw')
        data = GenerateGaussianData(sPreset.dim, sPreset.nGenDataCompnts, sPreset.N, sPreset.sDatasetParams.mu, sPreset.sDatasetParams.sigma, sPreset.sDatasetParams.p);
        sDataset.sData.x = data(1:sPreset.n,:);
        sDataset.sData.xt = data;
    end
    if sPreset.dim <= 3
        [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    end
    %figure; plot(sDataset.sData.x,sDataset.sData.y,'o',sDataset.sData.xt,sDataset.sData.yt,'.');
elseif strcmp(sPreset.verticesPDF, 'Uniform')
    if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
        sDataset.sData.x = GenerateUniformData(sPreset.dim, sPreset.n, sPreset.sDatasetParams.xMin, sPreset.sDatasetParams.xMax);
        sDataset.sData.xt = GenerateUniformData(sPreset.dim, sPreset.N, sPreset.sDatasetParams.xMin, sPreset.sDatasetParams.xMax);
    elseif strcmp(sPreset.dataGenTechnique, 'OneDraw')
        data = GenerateUniformData(sPreset.dim, sPreset.N, sPreset.sDatasetParams.xMin, sPreset.sDatasetParams.xMax);
        sDataset.sData.x = data(1:sPreset.n,:);
        sDataset.sData.xt = data;
    end
    [sDataset.sData.y, sDataset.sData.yt] = GenerateSyntheticGraphSignal(sDataset.sData.x, sDataset.sData.xt);
    
elseif strcmp(sPreset.verticesPDF, 'Grid')
    if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
        sDataset.sData.x = GenerateGridData(sPreset.dim, sPreset.sDatasetParams.nx, sPreset.sDatasetParams.xMin, sPreset.sDatasetParams.xMax);
        sDataset.sData.xt = GenerateGridData(sPreset.dim, sPreset.sDatasetParams.Nx, sPreset.sDatasetParams.xMin, sPreset.sDatasetParams.xMax);
    elseif strcmp(sPreset.dataGenTechnique, 'OneDraw')
        data = GenerateGridData(sPreset.dim, sPreset.sDatasetParams.Nx, sPreset.sDatasetParams.xMin, sPreset.sDatasetParams.xMax);
        msgBoxMsg = [];
        if sPreset.N ~= length(data)
            sPreset.N = length(data);
            msgBoxMsg = [msgBoxMsg, 'updated to sPreset.N = ' num2str(sPreset.N)];
        end
        r = ceil(sPreset.N/sPreset.n);
        if sPreset.n ~= length(1:r:sPreset.N)
            sPreset.n = length(1:r:sPreset.N);
            msgBoxMsg = [msgBoxMsg, 'updated to sPreset.n = ' num2str(sPreset.n)];
        end
        
        if sPreset.dim == 1
            xTrainInd = 1:r:sPreset.N;
        elseif sPreset.dim == 2
            assert(sqrt(r) == floor(sqrt(r)));
            [col, row] = meshgrid(1:sqrt(r):sPreset.sDatasetParams.Nx(1),1:sqrt(r):sPreset.sDatasetParams.Nx(2));
            xTrainInd = sub2ind(fliplr(sPreset.sDatasetParams.Nx),row(:),col(:));
        else
            error('invalid sPreset.dim for this dataset')
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
elseif strcmp(sPreset.verticesPDF, 'BulgariBeacons')
    sDataset = LoadBulgariBeacons(sPlotParams, sPreset.NLat, sPreset.NLon, sPreset.N, sPreset.n, sPreset.nLabeled);
else
    error('unknown pdf')
end

% Override y and yt with eigenvectors
if b_interpEigenvecs || (sPlotParams.b_globalPlotEnable && sPlotParams.b_plotWeights)
    [W, ~, dist, D, Ln, ~] = CalcAdjacency(sDataset.sData.x, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotWeights
        PlotWeightsMatrix([], W, dist, D, Ln, sDataset.sData.x, sPreset.adjacencyType, sPreset.omega, sPreset.k);
    end
end
if b_interpEigenvecs
    [WRef, ~, distRef, DRef, LnRef] = CalcAdjacency(sDataset.sData.xt, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
    [V, VRef] = EigsByTypeWrapper(sPlotParams, sPreset, sDataset, W, D, Ln, WRef, DRef, LnRef);
    interpRatio = sPreset.N/sPreset.n;
    if isfield(sDataset.sData, 'ymasked')
        sDataset.sData = rmfield(sDataset.sData, 'ymasked');
    end
    sDataset.sData.y = (1/sqrt(interpRatio))*V;
    sDataset.sData.yt = VRef;
end


% sanity checks
assert(sPreset.n == size(sDataset.sData.x,1) && sPreset.N == size(sDataset.sData.xt,1));
if strcmp(sPreset.dataGenTechnique, 'TwoDraws')
    warning('sPreset.sPreset.dataGenTechnique is TwoDraws, should make sure its okay...')
    pause(0.5)
else
    assert(isequal(sDataset.sData.x, sDataset.sData.xt(1:sPreset.n,:)));
end

fprintf('Done.\n')
end
