function [V, adjLambda, matLambda] = EigsByType(W, D, M, matrixForEigs)

if strcmp(matrixForEigs, 'Adjacency')
    [V, Lambda] = eigs(W,M);
    matLambda = diag(Lambda);
    adjLambda = matLambda;
elseif strcmp(matrixForEigs, 'RandomWalkDsvila')
    d = diag(D);
    S = diag(d.^-0.5)*W*diag(d.^-0.5);
    [V, Lambda] = eigs(S, M);
    matLambda = diag(Lambda);
    V = diag(d.^-0.5)*V;
    V = V./norm(V(:,1));
    
    adjLambda = eigs(W,M);
elseif strcmp(matrixForEigs, 'RandomWalk')
    % Reminder: the right and left eigenvectors of a non-symmetric matrix are not equal
    %   [Psi, Lambda, Phi] = eig(Pa) ---> 
    %       Pa*psi = lambda*psi       (Single column right eigenvector)
    %       phi^T*Pa = lambda*phi^T   (Single row left eigenvector)
    d = diag(D);
    alpha = 0;
    if alpha == 0
        Pa = D^-1*W;
    else
        Ka = diag(d.^-alpha)*W*diag(d.^-alpha);
        da = sum(Ka,2);
        Da = diag(da);
        Pa = Da^-1*Ka;
    end
    
    isalmostequal(diag(sum(Pa,2)),eye(length(Pa)),1e-13) % make sure Pa is row stochastic
    
    % Take the right eigenvectors of Pa
    [Psi, Lambda] = eigs(Pa, M); 
    matLambda = diag(Lambda);
    V = Psi;
    adjLambda = eigs(W,M);
elseif strcmp(matrixForEigs, 'NormLap')
    d = diag(D);
    Wn = diag(d.^-0.5)*W*diag(d.^-0.5);
    I = eye(length(W));
    Ln = I - Wn;
    [V, Lambda] = eig(Ln);
    V = V(:,1:M);
    matLambda = diag(Lambda);
    matLambda = matLambda(1:M);
    adjLambda = eigs(W,M);
elseif strcmp(matrixForEigs, 'Laplacian')
    L = D - W;
    [V, Lambda] = eig(L);
    V = V(:,1:M);
    matLambda = diag(Lambda);
    matLambda = matLambda(1:M);
    adjLambda = eigs(W,M);
else
    error('invalid matrixType = %s', matrixForEigs);
end

assert(~any(isnan(adjLambda)), 'You have at least one NaN');

end

