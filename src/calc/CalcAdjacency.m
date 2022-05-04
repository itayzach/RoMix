function [W, tW, dist, D, Ln, tLn] = CalcAdjacency(xTrain, adjacencyType, sDistanceParams, omega, k, nnValue)
fprintf('Calculating adjacency matrix for %d points... ',length(xTrain))
t = tic;
if strcmp(adjacencyType, 'GaussianKernel')
    n = size(xTrain,1);
    dist = CalcDistance(xTrain, xTrain, sDistanceParams);
    W = exp(-dist.^2/(2*omega^2));
    epsilon = 2*omega^2;
    
    distWithoutDiag = dist;
    distWithoutDiag(1:n+1:end) = [];
    distWithoutDiag = reshape(distWithoutDiag,n-1,n);
    assert(~isempty(find(distWithoutDiag > 0,1)), 'Maybe you have duplications of points?')
    minDist = min(distWithoutDiag);
    lafonEps = (1/n)*sum(minDist.^2);
    lafonOmega = sqrt(lafonEps/2);
%     fprintf('Lafon eps = %2.2f --> Lafon omega = 2*eps^2 = %2.2f\n', lafonEps, lafonOmega)
%     fprintf('Your  eps = %2.2f --> Your  omega = 2*eps^2 = %2.2f\n', epsilon, omega)
%     if(epsilon < max(min(distWithoutDiag.^2)))
%         warning(['Your graph is not fully connected, your omega = %.2f is too small. ',...
%          'Consider a bigger omega, like omega > %.2f\n'], ...
%             omega, 2*max(min(distWithoutDiag.^2))^2)
%     end
        
elseif strcmp(adjacencyType, 'NearestNeighbor')
    dist = []; % irrelevant for W
    W = NearestNeighborsAdjacency(xTrain, k, nnValue);
end
tW = toc(t);
fprintf('Done (took %.2f sec).\n', tW)

if nargout >= 4 
    fprintf('Calculating normalized Laplacian for %d points... ',length(xTrain))
    t = tic;
    d = sum(W,2);
    D = diag(d);
    Ln = eye(n) - diag(d.^-0.5)*W*diag(d.^-0.5);
    Ln = (Ln + Ln.')/2;
    tLn = toc(t);
    fprintf('Done (took %.2f sec).\n', tLn)
    fprintf('Checking connectivity by using eig(Ln)... ')
    if n < 4000
        t = tic;
        Ln_lambda = eig(Ln);
        nGraphComponents = find(Ln_lambda > 1e-5,1) - 1;
        assert(~isempty(nGraphComponents))
        fprintf('Done (took %.2f sec). G has %d connected components (Ln_lambda(%d) = %.3f)... ', toc(t), nGraphComponents, nGraphComponents+1, Ln_lambda(nGraphComponents+1))
        if nGraphComponents > 1
            distWithoutDiagSorted = sort(distWithoutDiag);
            ssl_heuristic_omega = (1/3)*mean(distWithoutDiagSorted(10,:));
            fprintf('\n')
            warning(['Your graph has %d connected components! Your omega = %.2f is too small. ', ...
                'Consider a bigger omega, like omega = %.2f\n'], nGraphComponents, omega ,ssl_heuristic_omega)
        end 
    else
        warning('n > 4000, eig will take time. skipping... ')
    end
end
end