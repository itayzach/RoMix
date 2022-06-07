function sDataset = LoadDigitsDataset(sPlotParams,actualDataDist, sDatasetParams, N, n)
    nSets = ceil(N/1000);
    datasetInd = randperm(10, nSets);
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
        [~, ~, ~, Ln] = CalcAdjacency(Xdiv255,'GaussianKernel',sPreset.sDistanceParams, sPreset.omega);
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

    
    if isfield(sDatasetParams, 'b_runVAE') && sDatasetParams.b_runVAE
        xTrainT = Xdiv255perm(1:n,:); 
        xIntT = Xdiv255perm;
        xTrainT = reshape(xTrainT,[],28,28);
        xTrain = permute(pagetranspose(permute(xTrainT,[2 3 1])),[3 1 2]);
        xTestT = reshape(xIntT(n+1:end,:),[],28,28);
        xTest = permute(pagetranspose(permute(xTestT,[2 3 1])),[3 1 2]);
        [zTrain, zTest, vae] = LoadMnistLatent(sPlotParams,sDatasetParams, xTrain, xTest);
        sDataset.vae = vae;
        sDataset.sData.x = zTrain;
        sDataset.sData.xt = [zTrain; zTest];

        sDataset.sData.y = mSignalsPerm(1:n,:);
        sDataset.sData.ymasked = mSignalsPerm(1:n,:);
        if nLabeled > n
            sDataset.sData.ymasked(nLabeled+1:end,:) = 0;
        end
        sDataset.sData.yt = mSignalsPerm; 
    else
        sDataset.sData.x = Xdiv255ordered(1:n,:);
        sDataset.sData.xt = Xdiv255ordered;
        sDataset.sData.y = mSignalsPermOrdered(1:n,:);
        sDataset.sData.ymasked = mSignalsPermOrderedMasked;
        sDataset.sData.yt = mSignalsPermOrdered;
    end
end