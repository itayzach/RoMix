function [mPhiNumeric, vLambdaNumeric] = CalcNumericEigenvectors(nEigs, sKernelParams, mData)

nTotal = length(mData);
nComponents = sKernelParams.sDistParams.estNumComponents;
%% CalcAdjacency
A = CalcAdjacency(sKernelParams, mData);

%% EVD
[mPhiNumeric, mLambdaNumeric] = eigs(A, nEigs);
[vLambdaNumeric, idx] = sort(diag(mLambdaNumeric), 'descend');

% Normalize
% vLambdaNumeric = (1/nTotal) * vLambdaNumeric;
% mPhiNumeric = sqrt(nTotal)*mPhiNumeric(:,idx);

vLambdaNumeric = (nComponents/nTotal)*vLambdaNumeric;
mPhiNumeric = mPhiNumeric(:,idx);

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
% mPhi_L = (n)*mPhi_L(:,1:sSimParams.CalcEigenFuncsM);


end