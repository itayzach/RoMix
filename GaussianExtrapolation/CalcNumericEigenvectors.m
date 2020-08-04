function [mPhi_A, vLambda_A] = CalcNumericEigenvectors(sSimParams, sKernelParams, mData)

n = length(mData);

%% CalcAdjacency
A = CalcAdjacency(sKernelParams,mData);

%% EVD
[mPhi_A, vLambda_A] = eigs(A, sSimParams.PlotSpectM);
[vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
vLambda_A = (1/n) * vLambda_A;
mPhi_A = sqrt(n)*mPhi_A(:,idx);

%% Laplacian
% D = diag(A*ones(n,1));
% t = sSimParams.t;
% d = sSimParams.dim;
% L = (((4*pi*t)^(-(d+2)/2))/n) * (D - A);
% 
% 
% [mPhi_L, vLambda_L] = eig(L);
% % [vLambda_L, idx] = sort(diag(vLambda_L), 'descend');
% vLambda_L = diag(vLambda_L);
% vLambda_L = vLambda_L(1:sSimParams.PlotSpectM);
% mPhi_L = (n)*mPhi_L(:,1:sSimParams.PlotEigenFuncsM);


end