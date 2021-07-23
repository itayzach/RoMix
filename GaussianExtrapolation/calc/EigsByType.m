function [V, adjLambda, matLambda] = EigsByType(W, M, matrixForEigs)

if strcmp(matrixForEigs, 'Adjacency')
    [V, Lambda] = eigs(W,M);
    matLambda = diag(Lambda);
elseif strcmp(matrixForEigs, 'RandomWalk')
%     d = sum(W,2);
%     S = diag(d.^-0.5)*W*diag(d.^-0.5);
%     [V, Lambda] = eigs(S, M);
%     matLambda = diag(Lambda);
%     V = diag(d.^-0.5)*V;
%     V = V./norm(V(:,1));

    alpha = 1;
    d = sum(W,2);
    Ka = diag(d.^-alpha)*W*diag(d.^-alpha);
    da = sum(Ka,2);
    Da = diag(da);
    Pa = Da^-1*Ka;
    isalmostequal(diag(sum(Pa,2)),eye(length(Pa)),1e-14) % make sure Pa is row stochastic
    [Psi, Lambda, Phi] = eigs(Pa,M);
    matLambda = diag(Lambda);
    V = Psi;
    % Psi = Psi./norm(Psi(:,1));
elseif strcmp(matrixForEigs, 'NormLap')
    d = sum(W,2);
    Wn = diag(d.^-0.5)*W*diag(d.^-0.5);
    I = eye(length(W));
    Ln = I - Wn;
    [V, Lambda] = eig(Ln);
    V = V(:,1:M);
    matLambda = diag(Lambda);
    matLambda = matLambda(1:M);
elseif strcmp(matrixForEigs, 'Laplacian')
    d = sum(W,2);
    D = diag(d);
    L = D - W;
    [V, Lambda] = eig(L);
    V = V(:,1:M);
    matLambda = diag(Lambda);
    matLambda = matLambda(1:M);

else
    error('invalid matrixType = %s', matrixForEigs);
end

adjLambda = eigs(W,M);

assert(~any(isnan(adjLambda)), 'You have at least one NaN');

end

