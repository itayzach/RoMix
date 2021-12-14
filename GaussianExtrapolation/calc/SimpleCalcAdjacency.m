function [W, dist, D, DsqrtInv] = SimpleCalcAdjacency(xTrain, adjacencyType, distType, omega, k, nnValue)
if strcmp(adjacencyType, 'GaussianKernel')
    n = length(xTrain);
    dist = CalcDistance(xTrain, xTrain, distType);
    W = exp(-dist.^2/(2*omega^2));
    epsilon = 2*omega^2;
    
    distWithoutDiag = dist;
    distWithoutDiag(1:n+1:end) = [];
    distWithoutDiag = reshape(distWithoutDiag,n-1,n);
    minDist = min(distWithoutDiag);
    recEps = (1/n)*sum(minDist.^2);
    recOmega = 2*recEps^2;
    assert(epsilon > max(min(distWithoutDiag.^2)), ...
        strcat('Your graph is not fully connected, your omega = %.2f is too small. ',...
            'Consider a bigger omega, like omega > %.2f\n'), ...
            omega, 2*max(min(distWithoutDiag.^2))^2)
elseif strcmp(adjacencyType, 'NearestNeighbor')
    dist = []; % irrelevant for W
    W = NearestNeighborsAdjacency(xTrain, k, nnValue);
end

d = sum(W,2);
D = diag(d);
DsqrtInv = diag(1./sqrt(d));
end