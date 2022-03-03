function [V, adjLambda, matLambda] = CalcAnalyticLapEigsGrid(x, M, len)
% analytic eigenfunctions from Spectral Graph Theory, Spielman
% are given by:
%     v_m(j) = cos(pi*m*j/n - pi*m/2n) 
%     lambda_m = 2(1-cos(pi*m/n))
% The analogy for us is by:
%     dx = L/n
%     x_j = j*dx - dx/2 = j*L/n - L/2n
%     ==> (pi*m/L)*x_j = pi*m*j/n - pi*m/2n 
[n, dim] = size(x);

V = zeros(n,M);
adjLambda = zeros(M,1);
OneDim2MultiDimIndexMatrix = LoadOneDim2MultiDimIndexMatrix(M, dim);
for i = 1:M
    m = OneDim2MultiDimIndexMatrix(i,:);
    normFactorV = (m == 0)*sqrt(1/n^(1/dim)) + (m > 0)*sqrt(2/n^(1/dim));
    V(:,i) = prod(normFactorV.*cos(pi*m.*x./len),2);
    adjLambda(i) = sum((pi*m./len).^2);
end
[adjLambda, vMultindexToSingleIndexMap] = sort(adjLambda);
V = V(:,vMultindexToSingleIndexMap);
matLambda = adjLambda;
end