function [V, lambda] = EigsByType(W, M, matrixForEigs)

if strcmp(matrixForEigs, 'Adjacency')
    [V, Lambda] = eigs(W,M);
    lambda = diag(Lambda);
elseif strcmp(matrixForEigs, 'RandomWalk')
%     d = sum(W,2);
%     S = diag(d.^-0.5)*W*diag(d.^-0.5);
%     [V, Lambda] = eigs(S, M);
%     lambda = diag(Lambda);
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
    lambda = diag(Lambda);
    V = Psi;
    % Psi = Psi./norm(Psi(:,1));
elseif strcmp(matrixForEigs, 'NormLap')
    d = sum(W,2);
    Wn = diag(d.^-0.5)*W*diag(d.^-0.5);
    I = eye(length(W));
    Ln = I - Wn;
    [V, Lambda] = eig(Ln);
    V = V(:,1:M);
    lambda = diag(Lambda);
    lambda = lambda(1:M);
elseif strcmp(matrixForEigs, 'Laplacian')
    d = sum(W,2);
    D = diag(d);
    L = D - W;
    [V, Lambda] = eig(L);
    V = V(:,1:M);
    lambda = diag(Lambda);
    lambda = lambda(1:M);

else
    error('invalid matrixType = %s', matrixForEigs);
end
assert(~any(isnan(lambda)), 'You have at least one NaN');

end

