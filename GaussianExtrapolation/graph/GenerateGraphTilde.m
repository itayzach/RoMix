function [G_tilde, dist_tilde, sDatasetTilde, sKernelParamsTilde] = GenerateGraphTilde(v_tilde,nTildeComponents,omega,nEigs,b_normlizedLaplacian)
sDatasetTilde.sData.x = v_tilde;
sDatasetTilde.estDataDist = 'Gaussian';
sDatasetTilde.dim = size(v_tilde,2);
GMMRegVal = 0;

% Estimate distribution and get kernel parameters
sDistParamsTilde = EstimateDistributionParameters(sDatasetTilde, nTildeComponents, GMMRegVal);
sKernelParamsTilde = GetKernelParams(sDatasetTilde, sDistParamsTilde, omega);
[sKernelParamsTilde.vLambdaAnaytic, sKernelParamsTilde.vComponentIndex, sKernelParamsTilde.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParamsTilde, sDatasetTilde.dim, nTildeComponents);
% Adjacency (gaussian kernel)
[W_tilde, dist_tilde] = CalcAdjacency(sKernelParamsTilde, v_tilde);
W_tilde = sparse(W_tilde);


% Graph
G_tilde = gsp_graph(W_tilde, v_tilde);
if b_normlizedLaplacian
    G_tilde.lap_type = 'normalized';
    G_tilde = gsp_graph_default_parameters(G_tilde);
end
end

