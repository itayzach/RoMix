function [ mPhiNys, vLambdaNys ] = CalcNystromEigenvectors(sParams, sKernelParams, mData, nysRatio)

nTotal = length(mData);

%% Split data matrix into blocks
mLeftUpperBlock = mData(1:nysRatio*nTotal,:);
mRightUpperBlock = mData(nysRatio*nTotal+1:end,:);

%% CalcAdjacency
A = CalcAdjacency(sKernelParams,mLeftUpperBlock);

%% EVD
[mPhi, vLambdaNys] = eigs(A, sParams.PlotSpectM);
[vLambdaNys, idx] = sort(diag(vLambdaNys), 'descend');
mPhi = mPhi(:,idx);

%% Nystrom
B = CalcAdjacency(sKernelParams,mLeftUpperBlock,mRightUpperBlock);
mPhiExt = B.'*mPhi*diag(1./vLambdaNys);
mPhiNys = sqrt(nTotal)*[mPhi; mPhiExt];
vLambdaNys = (1/nTotal) * vLambdaNys;
end