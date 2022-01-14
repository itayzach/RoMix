function sDistParams = EstimateDistributionParameters(x, gmmNumComponents, gmmRegVal, gmmMaxIter)

sDistParams.estDataDist = 'Gaussian';
dim = size(x,2);
sDistParams.dim = dim;

fprintf('Running fitgmdist with %d components... ', gmmNumComponents)
options = statset('MaxIter',gmmMaxIter);
GMModel = fitgmdist(x, gmmNumComponents, 'RegularizationValue', gmmRegVal, 'Options', options);
assert(GMModel.Converged, 'GMM couldn''t converge...')
sDistParams.GMModel = GMModel;
fprintf('Done.\n')


fprintf('Calculating [fliplr(u), fliplr(sigma^2)] = eig(cov) for all %d components... ',gmmNumComponents)
sDistParams.estNumComponents = gmmNumComponents;
sDistParams.componentProportion = GMModel.ComponentProportion;
for c = 1:gmmNumComponents
    sDistParams.cov{c} = GMModel.Sigma(:,:,c);
    sDistParams.mu{c} = GMModel.mu(c,:);
    [sDistParams.u{c}, sDistParams.sigma_eigv{c}] = eig(sDistParams.cov{c});
    sDistParams.sigma{c} = diag(sqrt(sDistParams.sigma_eigv{c})).';
    assert(isreal(sDistParams.sigma{c}) && all(sDistParams.sigma{c} > eps))
% 1:
%     sDistParams.u{c} = fliplr(sDistParams.u{c});
%     sDistParams.sigma{c} = fliplr(sDistParams.sigma{c});
% 2:
    [~, uIndByCov] = sort(sDistParams.sigma{c},'descend');
    sDistParams.u{c} = sDistParams.u{c}(:,uIndByCov);
    sDistParams.sigma{c} = sDistParams.sigma{c}(uIndByCov);

    sDistParams.mu_1D{c} = sDistParams.mu{c}*sDistParams.u{c};
    isalmostequal(sDistParams.u{c}*diag(sDistParams.sigma{c}.^2)*sDistParams.u{c}.', sDistParams.cov{c}, 1e-10)
%     assert(isequal(sort(sDistParams.sigma{c}),sDistParams.sigma{c}))
end
[minSigma, minSigmaInd] = min(sDistParams.sigma{c});
fprintf('Done. min(sigma{1:%d}) = %.4f (ind = %d)\n', gmmNumComponents, minSigma, minSigmaInd)

% Caluclate the probability for each data point x
% sDistParams.vPr = zeros(length(x),1);
% for c = 1:gmmNumComponents
%     prop = sDistParams.componentProportion(c);
%     vDensity = prop*p(sDistParams, c, x);
%     vDiffs = ones(length(x),1);
%     for d = 1:sDistParams.dim
%         [ xTotalSorted, vSortedIdx ] = sort(x(:,d));
%         [~, vInvSortIdx] = sort(vSortedIdx);
%         vDiffsSorted = [0; diff(xTotalSorted)];
%         vDiffs = vDiffs .* vDiffsSorted(vInvSortIdx);
%     end
%     sDistParams.vPr = sDistParams.vPr + (vDensity .* vDiffs);
% end

end