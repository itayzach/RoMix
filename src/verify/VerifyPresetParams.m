function VerifyPresetParams(sPreset, clusterMethod)
assert(~sPreset.b_debugUseAnalytic || ...
    (sPreset.b_debugUseAnalytic && ismember(sPreset.verticesPDF, {'Gaussian', 'Grid', 'Uniform', 'SwissRoll'})))
assert(~strcmp(sPreset.adjacencyType,'NearestNeighbor') || ...
    strcmp(sPreset.adjacencyType,'NearestNeighbor') && strcmp(sPreset.verticesPDF,'Grid'))
%assert((strcmp(clusterMethod, 'GMM') && sPreset.n >= sPreset.dim) || (strcmp(clusterMethod, 'SC')))
assert(strcmp(clusterMethod, 'GMM'))

if isfield(sPreset.sDatasetParams, 'nLabeled')
    assert(sPreset.sDatasetParams.nLabeled == sPreset.n);
end