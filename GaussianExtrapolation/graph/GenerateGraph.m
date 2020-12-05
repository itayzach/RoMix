function [G, dist, sDataset, sKernelParams] = GenerateGraph(graphName, nComponents,estDataDist,omega, b_normlizedLaplacian, nEigs, N)
%% Generate graph
switch graphName
    case {'Uniform_2D', 'Gaussian_2D', 'Uniform_1D', 'Gaussian_1D' }
        splitOut = split(graphName, "_");
        dataDist = splitOut{1};
        dataDim = str2double(splitOut{2}(1));
        %N = 1000;
        sDataset = GenerateDataset(dataDist, dataDim, nComponents, N, 0);
        GMMRegVal = 0;
    case 'Two_moons'
        %N = 1000;
        dataDist = 'Two_moons';
        dataDim = 2;
        sDataset = GenerateDataset(dataDist, dataDim, nComponents, N, 0);
        GMMRegVal = 0;
    case 'bunny'
        G_tmp = gsp_bunny();
    case 'twodgrid'
        %N = 500;
        N = round(sqrt(N))^2; 
        G_tmp = gsp_2dgrid(sqrt(N),sqrt(N)); 
    case 'sensor'
        %N = 500;
        G_tmp = gsp_sensor(N);
    case 'minnesota'
        if exist('N', 'var')
            error('number of nodes cannot be changed for minnesota graph');
        end
        G_tmp = gsp_minnesota();
    case 'david_sensor'
        assert(~exist('N', 'var'), 'number of nodes cannot be changed for david_sensor graph');
        G_tmp = gsp_david_sensor_network();
    case 'swiss_roll'
        %N = 500;
        G_tmp = gsp_swiss_roll(N);
    case 'random_ring'
        %N = 500;
        G_tmp = gsp_random_ring(N);
    otherwise
        error([graphName ' is not supported']);       
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
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
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

