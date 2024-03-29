function [V, matLambda] = CalcAnalyticLapEigsGrid(x, M, len)
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
matLambda = zeros(M,1);
OneDim2MultiDimIndexMatrix = LoadOneDim2MultiDimIndexMatrix(M, dim);
for m = 1:M
    I = OneDim2MultiDimIndexMatrix(m,:);
    normFactorV = 1; %(m == 0)*sqrt(1/n^(1/dim)) + (m > 0)*sqrt(2/n^(1/dim));
    V(:,m) = prod(normFactorV.*cos(pi*I.*x./len),2);
    matLambda(m) = sum((pi*I./len).^2);
end
[matLambda, vMultindexToSingleIndexMap] = sort(matLambda);
V = V(:,vMultindexToSingleIndexMap);
end