function [zTrain, zTest, vae] = LoadMnistLatent(sPlotParams,sDatasetParams,xTrain,xTest,b_plotDecoded)

if ~exist('b_plotDecoded', 'var')
    b_plotDecoded = true;
end

x_train = py.numpy.array(xTrain);
x_test = py.numpy.array(xTest);

%% Train / load trained model
if sDatasetParams.b_loadKeras
    modelName = 'keras_';
else
    modelName = [];
end
modelName = [modelName, 'd' num2str(sDatasetParams.latentDim), '_e', num2str(sDatasetParams.epochs), '_b', num2str(sDatasetParams.batchSize)];
if sDatasetParams.b_forceLoadTrainedVAE
    modelName = [modelName, '_n60000'];
else
    modelName = [modelName, '_n' num2str(size(xTrain,1))];
end
decSaveFolder = fullfile(pwd, 'vae','models',modelName,'dec_trained');
encSaveFolder = fullfile(pwd, 'vae','models',modelName,'enc_trained');
if isfolder(decSaveFolder) && isfolder(encSaveFolder)
    % Load from trained model
    vae = pyrunfile(fullfile("vae", "vae_load.py"), "vae", enc_save_folder=encSaveFolder, dec_save_folder=decSaveFolder);
    fprintf('Loaded VAE from files\n')
else
    % Train
    assert(~sDatasetParams.b_forceLoadTrainedVAE, 'You wanted to load %s, but you''re training...', modelName)
    MyMsgBox('Train VAE?', 'Train VAE', true);
    tic;
    [vae, history] = pyrunfile(fullfile("vae", "vae_train.py"), ["vae", "history"], x_train=x_train, x_test=x_test, ...
        latent_dim=uint32(sDatasetParams.latentDim), epochs=uint32(sDatasetParams.epochs), batch_size=uint32(sDatasetParams.batchSize), ...
        enc_save_folder=encSaveFolder, dec_save_folder=decSaveFolder);
    histMat = struct(history.history);
    figure; plot(cellfun(@double,cell(history.epoch)), cellfun(@double,cell(histMat.loss)));
    t = toc;
    fprintf('Training took %.2f minutes\n', t/60)
end


%% Encoder
[~, ~, z_train] = pyrunfile(fullfile("vae", "vae_encoder.py"), ["z_mean", "z_log_var", "z"], data=x_train,vae=vae);
[~, ~, z_test] = pyrunfile(fullfile("vae", "vae_encoder.py"), ["z_mean", "z_log_var", "z"], data=x_test,vae=vae);
zTrain = double(z_train);
zTest = double(z_test);

%% Plot decoded
if b_plotDecoded
    x_train_rec = pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=zTrain,vae=vae);
    xTrainRec = double(x_train_rec);
    plotInd = 1:50; %randperm(size(xTrain,1),50);
    xTrainRec = reshape(xTrainRec(plotInd,:,:),[],28*28);
    xTrainPlt = reshape(xTrain(plotInd,:,:),[],28*28);
    PlotDigits([],xTrainPlt,[],0,'Train','Train');
    PlotDigits(sPlotParams,xTrainRec,[],0,'Reconstructed','Reconstructed');
end

fprintf('Done.\n')
end