function A = NearestNeighborsAdjacency(X, k, nnValue)

[idx, D] = knnsearch(X, X, 'K', k);
N = length(X);
A = zeros(N,N);
for i = 1:N
    if strcmp(nnValue, 'ZeroOne')
        A(i,idx(i,:)) = 1;
        A(idx(i,:),i) = 1;
    elseif strcmp(nnValue, 'Distance')
        A(i,idx(i,:)) = D(i,:);
        A(idx(i,:),i) = D(i,:);
    else
        error('invalid nnValue')
    end
end
end