function [G, dist, sDataset, sKernelParams] = GenerateGraph(graphName, dataDim, nComponents,estDataDist,omega, b_normlizedLaplacian, nEigs)
%% Generate graph
switch graphName
    case {'Two_moons', 'Uniform', 'Gaussian'}
        N = 1000;
        sDataset = GenerateDataset(graphName, dataDim, nComponents, N, 0);
        GMMRegVal = 0;
    case 'bunny'
        G_tmp = gsp_bunny();
    case 'twodgrid'
        N = 500;
        N = round(sqrt(N))^2; G_tmp = gsp_2dgrid(sqrt(N),sqrt(N)); 
    case 'sensor'
        N = 500;
        G_tmp = gsp_sensor(N);
    case 'minnesota'
        G_tmp = gsp_minnesota();
    case 'david_seonsor'
        G_tmp = gsp_david_sensor_network();
    case 'swiss_roll'
        N = 500;
        G_tmp = gsp_swiss_roll(N);
    case 'random_ring'
        N = 500;
        G_tmp = gsp_random_ring(N);
    otherwise
        error([graphName 'is not supported']);       
end
    
if exist('G_tmp', 'var')
    sDataset.sData.x = G_tmp.coords;
    sDataset.estDataDist = estDataDist;
    sDataset.dim = size(sDataset.sData.x,2);
    GMMRegVal = 0.1;
end

%% Estimate distribution and get kernel parameters
sDistParams = EstimateDistributionParameters(sDataset, nComponents, GMMRegVal);
sKernelParams = GetKernelParams(sDataset, sDistParams, omega);
[sKernelParams.vLambdaAnaytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, sDataset.dim, nComponents);
v = sDataset.sData.x;
%% Adjacency (gaussian kernel)
[W, dist] = CalcAdjacency(sKernelParams, v);
W = sparse(W);
%% Graph
G = gsp_graph(W, v);
if b_normlizedLaplacian
    G.lap_type = 'normalized';
    G = gsp_graph_default_parameters(G);
end

end

