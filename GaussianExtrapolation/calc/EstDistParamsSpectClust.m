function sDistParams = EstDistParamsSpectClust(X, W, nClusters)

sDistParams.estDataDist = 'Gaussian';
dim = size(X,2);
sDistParams.dim = dim;

fprintf('Running spectralcluster...')
[clusterInd,V_SC,D_SC] = spectralcluster(W,nClusters,"Distance","precomputed");
fprintf(' done.\n')

% figure; scatter3(X(:,1), X(:,2), X(:,3),[],clusterInd, 'filled'); colormap(jet(nClusters)); colorbar;

sDistParams.GMModel.NumComponents = nClusters;
sDistParams.GMModel.SCcompIdx = clusterInd;
% figure; 
% for c = 1:nClusters
%     currentCompInd = clusterInd==c;
%     scatter3(X(currentCompInd,1), X(currentCompInd,2), X(currentCompInd,3),[],c*ones(1,sum(currentCompInd)), 'filled');
%     hold on;
% end

sDistParams.estNumComponents = nClusters;
for c = 1:nClusters
    X_c = X(clusterInd==c,:);
    nc = length(X_c);
    mu_c = mean(X_c); %sum(X_c,1)/nc; %mean(X_c);
    sDistParams.cov{c} = cov(X_c); %(nc-1)*(X_c - mu_c).'*(X_c - mu_c); %cov(X_c);
    sDistParams.mu{c} = mu_c;
    [sDistParams.u{c}, sDistParams.sigma_eigv{c}] = eig(sDistParams.cov{c});
    sDistParams.sigma{c} = diag(sqrt(sDistParams.sigma_eigv{c})).';
    if dim == 2
        if sDistParams.cov{c}(1,1) > sDistParams.cov{c}(2,2)
            warning('The variance in the first axis is greater than the variance in the second, but eig returns the eigenvalues in increasing order. So we fliplr')
            sDistParams.u{c} = fliplr(sDistParams.u{c});    
            sDistParams.sigma{c} = fliplr(sDistParams.sigma{c});
        end   
        sDistParams.u{c} = [-sDistParams.u{c}(:,1) -sDistParams.u{c}(:,2)];    
    end
    sDistParams.mu_1D{c} = sDistParams.mu{c}*sDistParams.u{c};
    isalmostequal(sDistParams.u{c}*diag(sDistParams.sigma{c}.^2)*sDistParams.u{c}.', sDistParams.cov{c}, 1e-10)
end


end