function [sDistParams, t] = EstimateDistributionParameters(x, gmmNumComponents, gmmRegVal, gmmMaxIter)

sDistParams.estDataDist = 'Gaussian';
dim = size(x,2);
sDistParams.dim = dim;

totalAttempts = 10;
fprintf('Running fitgmdist with %d components (%d attempts)... ', gmmNumComponents, totalAttempts)
converged = false;
numAttempts = 0;
warning('off','stats:gmdistribution:FailedToConverge');
while ~converged && numAttempts < totalAttempts
    ts = tic;
    options = statset('MaxIter',gmmMaxIter);
    GMModel = fitgmdist(x, gmmNumComponents, 'RegularizationValue', gmmRegVal, 'Options', options);
    t(1) = toc(ts);
    converged = GMModel.Converged;
    numAttempts = numAttempts + 1;
    fprintf('%d... ', numAttempts)
end
assert(converged, 'GMM couldn''t converge...')
fprintf('\nDone after %d attempts (%d iterations with AIC = %.2f, BIC =%.2f)\n', numAttempts, GMModel.NumIterations, GMModel.AIC, GMModel.BIC)

gmmNumComponents = GMModel.NumComponents;
sDistParams.GMModelPreOrder = GMModel;
fprintf('Calculating [u, sigma^2] = eig(cov) for all %d components (and sorting by sigma)... ',gmmNumComponents)
sDistParams.estNumComponents = gmmNumComponents;
ts = tic;
for c = 1:gmmNumComponents
    sDistParams.cov{c} = GMModel.Sigma(:,:,c);
    sDistParams.mu{c} = GMModel.mu(c,:);
    [sDistParams.u{c}, sDistParams.sigma_eigv{c}] = eig(sDistParams.cov{c});
    sDistParams.sigma{c} = diag(sqrt(sDistParams.sigma_eigv{c})).';
    assert(isreal(sDistParams.sigma{c}) && all(sDistParams.sigma{c} > eps))
    [~, uIndBySigma] = sort(sDistParams.sigma{c},'descend');
    sDistParams.u{c} = sDistParams.u{c}(:,uIndBySigma);
    sDistParams.sigma{c} = sDistParams.sigma{c}(uIndBySigma);

    sDistParams.mu_1D{c} = sDistParams.mu{c}*sDistParams.u{c};
    isalmostequal(sDistParams.u{c}*diag(sDistParams.sigma{c}.^2)*sDistParams.u{c}.', sDistParams.cov{c}, 1e-10)
%     assert(isequal(sort(sDistParams.sigma{c}),sDistParams.sigma{c}))
end
t(2) = toc(ts);

%sDistParams = SortComponentsByEigs(sDistParams, GMModel);
sDistParams.GMModel = GMModel;

[minSigmaAllComp, minSigmaCompDim] = cellfun(@min, sDistParams.sigma);
[minSigma, minSigmaCompInd] = min(minSigmaAllComp);
fprintf('Done (took %2.f sec). min(sigma{1:%d}) = %.4f (c = %d, dim = %d)\n', sum(t), gmmNumComponents, minSigma, minSigmaCompInd, minSigmaCompDim(minSigmaCompInd))

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

