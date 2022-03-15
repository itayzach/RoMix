function sDistParams = EstDistParamsSpectClust(X, W, nClusters)

sDistParams.estDataDist = 'Gaussian';
[n, dim] = size(X);
sDistParams.dim = dim;

fprintf('Running spectralcluster... ')
if dim > 3
    fprintf('Performing PCA... ')
    keyboard;
    [coeff,score,latent,tsquared,explained,mu] = pca(X);
    lowDim = explained > 0.05;
    XRec = score(:,lowDim)*coeff(:,lowDim)' + mu; %repmat(mu,13,1);
    [clusterInd, V_SC, D_SC] = spectralcluster(XRec,nClusters);
else
    %[clusterInd, V_SC, D_SC] = spectralcluster(W,nClusters,"Distance","precomputed");
    [clusterInd, V_SC, D_SC] = spectralcluster(X,nClusters);
end
fprintf('Done.\n')

sDistParams.SCcompIdx = clusterInd;

if dim == 3
    fig = figure; ax = scatter3(X(:,1), X(:,2), X(:,3),[],clusterInd, 'filled'); colormap(jet(nClusters)); colorbar;
    UpdateCursorDataTip(fig, ax, clusterInd);
elseif dim == 2
    figure; scatter3(X(:,1), X(:,2), clusterInd,[],clusterInd, 'filled'); colormap(jet(nClusters)); colorbar; view(2)
end
if nClusters > 1
    figure; 
    scatter3(V_SC(:,2), V_SC(:,3), clusterInd,[],clusterInd, 'filled'); 
    %view(45,45)
    colormap(jet(nClusters))
    h = colorbar('TickLabelInterpreter', 'latex');
    h.TickLabels = [1 dim];
    h.Ticks = (1+0.25):0.5:(clusterInd-0.25);
    figure; histogram(clusterInd)
end
% for c = 1:nClusters
%     X_c = X(clusterInd==c,:);
%     figure; 
%     tiledlayout('flow')
%     imsize = sqrt(size(X_c,2));
%     for i = 1:size(X_c,1)
%         nexttile;
%         image = reshape(X_c(i,:),imsize,imsize) ;
%         imagesc(image)
%     end
%     sgtitle('random samples')
% end

fprintf('Calculating covariance, mean and [fliplr(u), fliplr(sigma^2)] = eig(cov) for all %d components... ',nClusters)
sDistParams.estNumComponents = nClusters;
for c = 1:nClusters
    X_c = X(clusterInd==c,:);
    assert(size(X_c,1) > 1, 'cluster %d has only one point',c)
    mu_c = mean(X_c,1);
    cov_c = cov(X_c);
    assert(det(cov_c) > 0, 'det(cov_%d) < 0?', c);
    [u_c, sigma_eigv_c] = eig(cov_c);
    assert(isreal(sigma_eigv_c), 'sigma_eigv_%d is complex?', c)
    assert(isreal(u_c), 'u_%d is complex?', c)
    assert(sum(sigma_eigv_c(:)) > 0, 'all sigma_eigv_c are zero!')
    
    sDistParams.cov{c} = cov_c;
    sDistParams.mu{c} = mu_c;
    %sDistParams.GMModel.mu(c,:) = mu_c;
    sDistParams.u{c} = u_c;
    sDistParams.sigma_eigv{c} = sigma_eigv_c;
    sDistParams.sigma{c} = diag(sqrt(abs(sigma_eigv_c))).';
    if ~all(sDistParams.sigma{c} > eps)
        warning('some of the sigmas are zero. sum(sigma>eps) = %d', sum(sDistParams.sigma{c} > eps))
    end
    assert(sum(sDistParams.sigma{c} > 0), 'all sigmas are zero!')

    [~, uIndByCov] = sort(sDistParams.sigma{c},'descend');
    sDistParams.u{c} = sDistParams.u{c}(:,uIndByCov);
    sDistParams.sigma{c} = sDistParams.sigma{c}(uIndByCov);
    sDistParams.mu_1D{c} = sDistParams.mu{c}*sDistParams.u{c};
    
    isalmostequal(sDistParams.u{c}*diag(sDistParams.sigma{c}.^2)*sDistParams.u{c}.', sDistParams.cov{c}, 1e-10)
end

mMu = cell2mat(sDistParams.mu');
mCov = reshape(cell2mat(sDistParams.cov), dim, dim, nClusters);
p = histcounts(clusterInd, nClusters)/n;
GMModel = gmdistribution(mMu, mCov, p);
sDistParams = SortComponentsByEigs(sDistParams, GMModel);

[minSigmaAllComp, minSigmaCompDim] = cellfun(@min, sDistParams.sigma);
[minSigma, minSigmaCompInd] = min(minSigmaAllComp);

fprintf('Done. min(sigma{1:%d}) = %.4f (c = %d, dim = %d)\n', nClusters, minSigma, minSigmaCompInd, minSigmaCompDim(minSigmaCompInd))

end