function sDistParams = SortComponentsByEigs(sDistParams, GMModel)
mCovEigs = cell2mat(sDistParams.sigma');
minEig = min(mCovEigs(:));
% Sum how many eigenvaluess are more than 5% greater than the minimum eigenvalue,
% and sort the components by that number
vHowManyMoreThanMinComp = sum(mCovEigs > 1.05*minEig,2);
%vHowManyMoreThanMinComp = sum(mCovEigs,2);
[~, vPrincipalCompsOrder] = sort(vHowManyMoreThanMinComp,'descend');


% Sort
sDistParams.cov = sDistParams.cov(vPrincipalCompsOrder);
sDistParams.mu = sDistParams.mu(vPrincipalCompsOrder);
sDistParams.u = sDistParams.u(vPrincipalCompsOrder);
sDistParams.sigma_eigv = sDistParams.sigma_eigv(vPrincipalCompsOrder);
sDistParams.sigma = sDistParams.sigma(vPrincipalCompsOrder);
sDistParams.mu_1D = sDistParams.mu_1D(vPrincipalCompsOrder);
sDistParams.GMModel = gmdistribution(GMModel.mu(vPrincipalCompsOrder,:), ...
    GMModel.Sigma(:,:,vPrincipalCompsOrder), GMModel.ComponentProportion(vPrincipalCompsOrder));
end