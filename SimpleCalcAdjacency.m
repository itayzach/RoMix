function W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega, k)
if strcmp(origGraphAdjacency, 'GaussianKernel')
    dist = pdist2(xTrain, xTrain);
    W = exp(-dist.^2/(2*omega^2));
elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
    error('Nystrom cant work with this type')
    W = NearestNeighborsAdjacency(xTrain, k);
end