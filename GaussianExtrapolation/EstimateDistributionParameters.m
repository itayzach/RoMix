function sDistParams = EstimateDistributionParameters(sDataset)

sDistParams.estDataDist = sDataset.estDataDist;
sDistParams.dim = sDataset.dim;

if strcmp(sDataset.estDataDist, 'Gaussian')
    GMModel = fitgmdist(sDataset.sData.x,1);
    sDistParams.cov = GMModel.Sigma;
    sDistParams.mu  = GMModel.mu;

    [sDistParams.u, sDistParams.sigma_eigv] = eig(sDistParams.cov);
    sDistParams.sigma = diag(sqrt(sDistParams.sigma_eigv)).';
    if sDataset.dim == 2
        if sDistParams.cov(1,1) > sDistParams.cov(2,2)
            warning('The variance in the first axis is greater than the variance in the second, but eig returns the eigenvalues in increasing order. So we fliplr')
            sDistParams.u = fliplr(sDistParams.u);    
            sDistParams.sigma = fliplr(sDistParams.sigma);
        end   
        sDistParams.u = [-sDistParams.u(:,1) -sDistParams.u(:,2)];    
    end
    sDistParams.mu_1D = sDistParams.mu*sDistParams.u;
    isalmostequal(sDistParams.u*diag(sDistParams.sigma.^2)*sDistParams.u.', sDistParams.cov, 1e-15)

    sDistParams.xMax = sDistParams.mu + 1.5*sDistParams.sigma;
    sDistParams.xMin = sDistParams.mu - 1.5*sDistParams.sigma;
    if sDistParams.xMax - sDistParams.xMin > 10
        sDistParams.xMax = sDistParams.mu + 5;
        sDistParams.xMin = sDistParams.mu - 5;
        warning('3sigma is too large...');
    end
    
    % Caluclate the probability for each data point x
    if sDataset.dim == 1
        [ xTotalSorted, vSortedIdx ] = sort(sDataset.sData.x);
        [~, vInvSortIdx] = sort(vSortedIdx);
        vDiffsSorted = [0; diff(xTotalSorted)];
        vDiffs = vDiffsSorted(vInvSortIdx);
        vDensity = p(sDistParams, sDataset.sData.x, sDataset.dim);
        sDistParams.vPr = vDensity .* vDiffs;
    end
    
elseif strcmp(sDataset.estDataDist, 'Uniform')
    sDistParams.a    = 2.5;
    sDistParams.u    = eye(sDistParams.dim);
    sDistParams.xMax = sDistParams.a;
    sDistParams.xMin = -sDistParams.a;
else
    error('unknown pdf')
end
end