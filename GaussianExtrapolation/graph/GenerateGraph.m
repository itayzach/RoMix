function [graphName, G, dist, sDataset, sKernelParams] = GenerateGraph(sSimParams, b_distributionDataset,nComponents,estDataDist,omega, b_normlizedLaplacian, nEigs)
%% Generate graph
if b_distributionDataset
    N = 1000;
    dataDim = 1;
    graphName = 'Uniform';
    sDataset = GenerateDataset(graphName, dataDim, nComponents, N, 0);
    GMMRegVal = 0;
else
    N = 500;
%     N = round(sqrt(N))^2; G_tmp = gsp_2dgrid(sqrt(N),sqrt(N)); graphName = 'twodgrid';
%     G_tmp = gsp_sensor(N);  graphName = 'sensor';
%     G_tmp = gsp_bunny(); N = G_tmp.N; graphName = 'bunny';
    G_tmp = gsp_minnesota(); N = G_tmp.N; graphName = 'minnesota';
%     G_tmp = gsp_david_sensor_network(); N = G_tmp.N; graphName = 'david_seonsor';
%     G_tmp = gsp_swiss_roll(N); graphName = 'swiss_roll';
%     G_tmp = gsp_random_ring(N); graphName = 'random_ring';
    sDataset.sData.x = G_tmp.coords;
    sDataset.estDataDist = estDataDist;
    sDataset.dim = size(sDataset.sData.x,2);
    GMMRegVal = 0.1;
    dataDim = size(G_tmp.coords,2);
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

%% Plot the graph
if sSimParams.b_showGraphMatricesFigures
    figure; 
    subplot(131)
        imagesc(G.W); 
        colorbar; 
        title('Kerel (weights) matrix', 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
    subplot(132); 
        imagesc(diag(G.d)); 
        colorbar; 
        title('Degree matrix', 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);      
    subplot(133); 
        imagesc(G.L); 
        colorbar; 
        if b_normlizedLaplacian
            title('Normalized Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
        else
            title('Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
        end
        set(gca,'FontSize', 14);    
    set(gcf,'Position', [100 200 1800 400])
end
end

