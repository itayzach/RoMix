function [W, dist, D] = SimpleCalcAdjacency(xTrain, adjacencyType, omega, k, nnValue)
if strcmp(adjacencyType, 'GaussianKernel')
    dist = pdist2(xTrain, xTrain);
    W = exp(-dist.^2/(2*omega^2));
elseif strcmp(adjacencyType, 'NearestNeighbor')
    dist = []; % irrelevant for W
    W = NearestNeighborsAdjacency(xTrain, k, nnValue);
end
D = diag(sum(W,2));
end