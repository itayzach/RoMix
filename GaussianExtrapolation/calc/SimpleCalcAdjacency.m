function [W, dist, D, Ln] = SimpleCalcAdjacency(xTrain, adjacencyType, distType, omega, k, nnValue)
if strcmp(adjacencyType, 'GaussianKernel')
    n = length(xTrain);
    dist = CalcDistance(xTrain, xTrain, distType);
    W = exp(-dist.^2/(2*omega^2));
    epsilon = 2*omega^2;
    
    distWithoutDiag = dist;
    distWithoutDiag(1:n+1:end) = [];
    distWithoutDiag = reshape(distWithoutDiag,n-1,n);
    minDist = min(distWithoutDiag);
    lafonEps = (1/n)*sum(minDist.^2);
    lafonOmega = sqrt(lafonEps/2);
    fprintf('Lafon eps = %2.2f --> Lafon omega = 2*eps^2 = %2.2f\n', lafonEps, lafonOmega)
    fprintf('Your  eps = %2.2f --> Your  omega = 2*eps^2 = %2.2f\n', epsilon, omega)
    if(epsilon < max(min(distWithoutDiag.^2)))
        warning(['Your graph is not fully connected, your omega = %.2f is too small. ',...
         'Consider a bigger omega, like omega > %.2f\n'], ...
            omega, 2*max(min(distWithoutDiag.^2))^2)
    end
        
elseif strcmp(adjacencyType, 'NearestNeighbor')
    dist = []; % irrelevant for W
    W = NearestNeighborsAdjacency(xTrain, k, nnValue);
end

d = sum(W,2);
D = diag(d);
Ln = eye(n) - diag(d.^-0.5)*W*diag(d.^-0.5);
end