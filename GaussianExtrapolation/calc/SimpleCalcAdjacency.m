function [W, dist] = SimpleCalcAdjacency(xTrain, adjacencyType, omega, k, nnValue)
if strcmp(adjacencyType, 'GaussianKernel')
    dist = pdist2(xTrain, xTrain);
    W = exp(-dist.^2/(2*omega^2));
elseif strcmp(adjacencyType, 'NearestNeighbor')
%     error('Nystrom cant work with this type')
    dist = [];
    W = NearestNeighborsAdjacency(xTrain, k, nnValue);
end