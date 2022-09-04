clear; clc; close all;

tic;
mnistFromKeras = 0;
latentDim = 20;
epochs = 150;
batchSize = 256;
%nSamples = 5000;
%outputArgs = ["vae", "x_train","x_test"];
%[vae, x_train, x_test] = pyrunfile("vae_train.py", outputArgs, latent_dim=uint32(latentDim),num_samples=uint32(nSamples));

sPreset = GetMnistPreset();
nTrain = sPreset.n;
nTest = sPreset.N - sPreset.n;
if mnistFromKeras
    [x_train, x_test, y_train, y_test] = pyrunfile("load_mnist_keras.py", ["x_train", "x_test", "y_train", "y_test"], num_train=uint32(nTrain),num_test=uint32(nTest));
    xTrain = double(x_train);
    xTest = double(x_test);
    yTrain = double(y_train).';
    yTest = double(y_test).';
else
    sDataset = GenerateDataset(sPreset.verticesPDF, sPreset.dim, sPreset.nGenDataCompnts, sPreset.n, sPreset.N, sPreset.dataGenTechnique, sPreset.sDatasetParams);
    xTrainT = sDataset.sData.x; yTrain = sDataset.sData.y; xIntT = sDataset.sData.xt;
    xTrainT = reshape(xTrainT,[],28,28);
    xTrain = permute(pagetranspose(permute(xTrainT,[2 3 1])),[3 1 2]);
    xTestT = reshape(xIntT(sPreset.n+1:end,:),[],28,28);
    xTest = permute(pagetranspose(permute(xTestT,[2 3 1])),[3 1 2]);
    x_train = py.numpy.array(xTrain);
    x_test = py.numpy.array(xTest);
end
[vae, history] = pyrunfile("vae_train.py", ["vae", "history"], x_train=x_train, x_test=x_test, latent_dim=uint32(latentDim), epochs=uint32(epochs), batch_size=uint32(batchSize));
histMat = struct(history.history); 
figure; plot(cellfun(@double,cell(history.epoch)), cellfun(@double,cell(histMat.loss)));

d_x_train = double(x_train);
%x_test = double(x_test);

%% enc
[z_mean, z_log_var, z_train] = pyrunfile("vae_encoder.py", ["z_mean", "z_log_var", "z"], data=x_train,vae=vae);
zTrain = double(z_train);
PlotDataset([], [], zTrain, yTrain, 'plot');
x_train_rec = pyrunfile("vae_decoder.py", "x", z=zTrain,vae=vae);
xTrainRec = double(x_train_rec);
plotInd = randperm(nTrain,50);
xTrainRec = reshape(xTrainRec(plotInd,:,:),[],28*28);
PlotDigits([],reshape(xTrain(plotInd,:,:),[],28*28),[],0,'Train');
PlotDigits([],xTrainRec,[],0,'Reconstructed');
%% GMM
%nGmmPoints = 5000;
%zz = py.numpy.array(randn(nGmmPoints,latentDim));
%d_zz = double(zz);

gmmNumComponents = 25; gmmRegVal = 1e-5; gmmMaxIter = 2000;
sDistParams = EstimateDistributionParameters(zTrain, gmmNumComponents, gmmRegVal, gmmMaxIter);

% PlotGaussianEllipses([], sDistParams);
PlotCovEigs([], sDistParams);
nGmmPoints = nTrain;
[zGmm,compIdx] = random(sDistParams.GMModel, nGmmPoints);
PlotDataset([], [], zGmm, [], 'gmm');
zGmm = py.numpy.array(zGmm);
x = pyrunfile("vae_decoder.py", "x", z=zGmm,vae=vae);
d_x = double(x); 
d_x = reshape(d_x(1:50,:,:),[],28*28);
PlotDigits([],d_x,[],0,'GMM');
toc;



