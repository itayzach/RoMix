function [ mPhiNys, vLambdaNys ] = CalcNystromEigenvectors(sSimParams, sKernelParams, mData, nysRatio)

nTotal = length(mData);
nComponents = sKernelParams.sDistParams.estNumComponents;
%% Split data matrix into blocks
mLeftUpperBlock = mData(1:nysRatio*nTotal,:);
mRightUpperBlock = mData(nysRatio*nTotal+1:end,:);

%% CalcAdjacency
A = CalcAdjacency(sKernelParams, mLeftUpperBlock);
B = CalcAdjacency(sKernelParams, mLeftUpperBlock, mRightUpperBlock);

%% EVD with no orthogonalization and Nystrom
[mPhi, mLambdaNys] = eigs(A, sSimParams.PlotSpectM);
[vLambdaNys, idx] = sort(diag(mLambdaNys), 'descend');
mPhi = mPhi(:,idx);
mPhiExt = B.'*mPhi*diag(1./vLambdaNys);

% Normalize
% vLambdaNys = (1/nTotal) * vLambdaNys;
% mPhiNys = sqrt(nTotal)*[mPhi; mPhiExt];

vLambdaNys = (nComponents/nTotal)*vLambdaNys;
mPhiNys = [mPhi; mPhiExt];

%% EVD with orthgonalization (from Fowlkes)
% **************************************
% *** Pretty heavy computationaly... ***
% **************************************
% n = nysRatio*nTotal;
% m = round((1-nysRatio)*nTotal);
% d1 = sum([A;B.'],1);
% d2 = sum(B,1) + sum(B.',1)*pinv(A)*B;
% dhat = sqrt(1./[d1 d2]).';
% A = A.*(dhat(1:n)*dhat(1:n).');
% B = B.*(dhat(1:n)*dhat(n+(1:m)).');
% Asi = real(sqrtm(pinv(A)));
% S = A + Asi*(B*B.')*Asi;
% 
% % [U, mLambdaS] = eigs(S, sSimParams.PlotSpectM);
% % [vLambdaNys, idx] = sort(diag(mLambdaS), 'descend');
% % U = U(:,idx);
% 
% [U, mSigmaS, ~] = svd(S);
% U = U(:,1:sSimParams.CalcEigenFuncsM);
% vLambdaNys = sqrt(diag(mSigmaS));
% vLambdaNys = vLambdaNys(1:sSimParams.CalcEigenFuncsM);
% 
% mPhiNys = [A; B.']*Asi*U*(diag(1./vLambdaNys));
% 
% % mPhiNys = sqrt(nTotal)*mPhiNys;
% % vLambdaNys = (1/nTotal)*vLambdaNys;
