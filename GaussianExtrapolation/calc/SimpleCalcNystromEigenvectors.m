function [VNys, lambdaNys] = SimpleCalcNystromEigenvectors(x, omega, M, nysRatio)
%% Split data matrix into blocks
nTotal = length(x);
LUBlock = x(1:round(nysRatio*nTotal),:);
RUBlock = x(round(nysRatio*nTotal)+1:end,:);

%% CalcAdjacency - A
distLUBlock = pdist2(LUBlock, LUBlock);
A = exp(-distLUBlock.^2/(2*omega^2));

%% CalcAdjacency - B
% distLUBlockRUBlock = pdist2(LUBlock, RUBlock);
distLUBlockRUBlock = pdist2(LUBlock, x);
B = exp(-distLUBlockRUBlock.^2/(2*omega^2));

%% EVD with no orthogonalization and Nystrom
[V, LambdaNys] = eigs(A, M);
[lambdaNys, idx] = sort(diag(LambdaNys), 'descend');
V = V(:,idx);
VNys = B.'*V*diag(1./lambdaNys);
end
