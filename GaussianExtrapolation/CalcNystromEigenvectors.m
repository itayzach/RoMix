function [ mPhiNys, vLambdaNys ] = CalcNystromEigenvectors(sParams)

x = sParams.x_rand(1:sParams.nPointsNystrom,:);
y = sParams.x_rand(sParams.nPointsNystrom+1:end,:);
n = length(sParams.x_rand);

%% CalcAdjacency
A = CalcAdjacency(sParams,x);

%% EVD
[mPhi, vLambdaNys] = eigs(A, sParams.PlotSpectM);
[vLambdaNys, idx] = sort(diag(vLambdaNys), 'descend');
mPhi = mPhi(:,idx);

%% Nystrom
B = CalcAdjacency(sParams,x,y);
mPhiExt = B.'*mPhi*diag(1./vLambdaNys);
mPhiNys = sqrt(n)*[mPhi; mPhiExt];
vLambdaNys = (1/n) * vLambdaNys;
end