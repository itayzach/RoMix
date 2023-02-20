function sDataset = LoadDigitsKeras(sPlotParams,sDatasetParams, N, n, nLabeled)
nTrain = n;
nTest = N - n;
assert(nTest <= 10000, 'nTest should be less than 10000')
[x_train, x_test, y_train, y_test] = pyrunfile(fullfile("vae", "load_mnist_keras.py"), ["x_train", "x_test", "y_train", "y_test"], num_train=uint32(nTrain),num_test=uint32(nTest));
xTrain = double(x_train);
xTest = double(x_test);
yTrain = double(y_train).';
yTest = double(y_test).';

oneHotTrain = false(nTrain,10);
for i = 1:10
    oneHotTrain(yTrain==i-1,i) = true;
end
oneHotTest = false(nTest,10);
for i = 1:10
    oneHotTest(yTest==i-1,i) = true;
end

% Zero out labeled signals for ssl
vNonLabeledInd = randperm(nTrain,nTrain-nLabeled);

if isfield(sDatasetParams, 'b_runVAE') && sDatasetParams.b_runVAE
    [zTrain, zTest, vae] = LoadMnistLatent(sPlotParams,sDatasetParams, xTrain, xTest);
    sDataset.vae = vae;
    sDataset.sData.x = zTrain;
    sDataset.sData.xt = [zTrain; zTest];

    sDataset.sData.y = oneHotTrain;
    sDataset.sData.ymasked = oneHotTrain;
    sDataset.sData.ymasked(vNonLabeledInd,:) = 0;
    sDataset.sData.yt = [oneHotTrain; oneHotTest];
else
    xTrain = reshape(permute(xTrain, [1 3 2]),nTrain,[]);
    xTest =  reshape(permute(xTest, [1 3 2]),nTest,[]);
    sDataset.sData.x = xTrain;
    sDataset.sData.xt = [xTrain; xTest];
    sDataset.sData.y = oneHotTrain;
    sDataset.sData.yt = [oneHotTrain; oneHotTest];
    sDataset.sData.ymasked = oneHotTrain;
    sDataset.sData.ymasked(vNonLabeledInd,:) = 0;
end

end