function sDistParams = EstimateDistributionParameters(x, estNumComponents, gmmRegVal)

sDistParams.estDataDist = 'Gaussian';
dim = size(x,2);
sDistParams.dim = dim;

lastwarn('')
try
    GMModel = fitgmdist(x, estNumComponents, 'RegularizationValue', gmmRegVal);
%     GMModel = fitgmdist(x, estNumComponents, 'RegularizationValue', 0, 'SharedCovariance', true);
catch exception
    sprintf('GM error: %s\n', exception.message)
end
[warnMsg, ~] = lastwarn;
if ~isempty(warnMsg)
    sprintf('GM error: %s\n', warnMsg)
end
assert(GMModel.Converged)
sDistParams.GMModel = GMModel;

sDistParams.estNumComponents = estNumComponents;
sDistParams.componentProportion = GMModel.ComponentProportion;
for c = 1:estNumComponents
    sDistParams.cov{c} = GMModel.Sigma(:,:,c);
    sDistParams.mu{c} = GMModel.mu(c,:);
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

% Caluclate the probability for each data point x
sDistParams.vPr = zeros(length(x),1);
for c = 1:estNumComponents
    prop = sDistParams.componentProportion(c);
    vDensity = prop*p(sDistParams, c, x);
    vDiffs = ones(length(x),1);
    for d = 1:sDistParams.dim
        [ xTotalSorted, vSortedIdx ] = sort(x(:,d));
        [~, vInvSortIdx] = sort(vSortedIdx);
        vDiffsSorted = [0; diff(xTotalSorted)];
        vDiffs = vDiffs .* vDiffsSorted(vInvSortIdx);
    end
    sDistParams.vPr = sDistParams.vPr + (vDensity .* vDiffs);
end

end