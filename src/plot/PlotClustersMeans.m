function PlotClustersMeans(sPreset, sPlotParams, sDistParams)
if ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})
    vSamples = 1:sDistParams.estNumComponents;
    if isfield(sDistParams, 'vae')
        % Transform from z to x
        b_transpose = false;
        x = cell2mat(sDistParams.mu');
        x_train_rec = pyrunfile(fullfile("vae", "vae_decoder.py"), "x", z=x,vae=sDistParams.vae);
        xTrainRec = double(x_train_rec);
        xMeans = reshape(xTrainRec,[],28*28);
    else
        xMeans = cell2mat(sDistParams.mu');
        b_transpose = strcmp(sPreset.verticesPDF, 'MNIST');
    end
    
    PlotDigits(sPlotParams, xMeans(vSamples,:), vSamples, b_transpose, 'Means', 'Means');
end
% mMu = cell2mat(sDistParams.mu');
% mSigma = cell2mat(sDistParams.sigma');
% nComp = sDistParams.estNumComponents;
% dim = sDistParams.dim;
% mU = reshape(cell2mat(sDistParams.u),dim,dim,nComp);
end