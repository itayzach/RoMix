function vPr = EstimateProbabilityArea(sDistParams, x)

vPr = zeros(length(x),1);
for c = 1:sDistParams.estNumComponents
    prop = sDistParams.componentProportion(c);
    vDensity = prop*p(sDistParams, c, x);
    vDiffs = ones(length(x),1);
    for d = 1:sDistParams.dim
        [ xTotalSorted, vSortedIdx ] = sort(x(:,d));
        [~, vInvSortIdx] = sort(vSortedIdx);
        vDiffsSorted = [0; diff(xTotalSorted)];
        vDiffs = vDiffs .* vDiffsSorted(vInvSortIdx);
    end
    vPr = vPr + (vDensity .* vDiffs);
end
end