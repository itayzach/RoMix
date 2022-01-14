function OneDim2MultiDimIndexMatrix = LoadOneDim2MultiDimIndexMatrix(nEigs,dim)
fname = ['OneDim2MultiDimIndexMatrix_',num2str(dim),'d.mat'];
fpath = fullfile('GaussianExtrapolation','indmat');
if ~isfolder(fpath)
    fprintf('Creating %s... ',fpath)
    mkdir(fpath)
end
if isfile(fullfile(fpath,fname))
    fprintf('Loading %s...', fname)
    load(fullfile(fpath,fname), 'OneDim2MultiDimIndexMatrix')
end
if ~isfile(fullfile(fpath,fname)) || (nEigs > size(OneDim2MultiDimIndexMatrix,1))
    fprintf('\n')
    warning('Generating %s (this might take some time)... ', fname)
    OneDim2MultiDimIndexMatrix = zeros(nEigs,dim);
    for i=1:nEigs
        OneDim2MultiDimIndexMatrix(i,:) = OneDim2MultiDimIndex(i-1,dim);
    end
    save(fullfile(fpath,fname));
    fprintf('\n')
end
end