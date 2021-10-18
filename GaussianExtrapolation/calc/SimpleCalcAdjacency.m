function [W, dist, D] = SimpleCalcAdjacency(xTrain, adjacencyType, distType, omega, k, nnValue)
if strcmp(adjacencyType, 'GaussianKernel')
    dist = CalcDistance(xTrain, xTrain, distType);
    W = exp(-dist.^2/(2*omega^2));
%     W = W - eye(length(W));
%     W(W<0.999) = 0;
%     W = W/W(1,2);
    epsilon = 2*omega^2;
    distForAssert = dist;
    distForAssert(1:length(xTrain):end) = inf;
    assert(epsilon > max(min(distForAssert.^2)), ...
        'Your graph is not connected. Your omega = %.2f. Consider omega > %.2f\n', ...
        omega, 2*max(min(distForAssert.^2))^2)
elseif strcmp(adjacencyType, 'NearestNeighbor')
    dist = []; % irrelevant for W
    W = NearestNeighborsAdjacency(xTrain, k, nnValue);
%     W(1,k) = 0;
%     W(k,1) = 0;
%     W(end,end-k+1) = 0;
%     W(end-k+1,end) = 0;
end

D = diag(sum(W,2));
end