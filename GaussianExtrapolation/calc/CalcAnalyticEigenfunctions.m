function [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, mData, b_normalize)

[nTotal, dim] = size(mData);
fprintf('Calculating %d eigenfunctions (b_normalize = %d) d = %d, n = %d... ',nEigs,b_normalize, dim,nTotal)
t = tic;
mPhiAnalytic = zeros(nTotal, nEigs);
OneDim2MultiDimIndexMatrix = LoadOneDim2MultiDimIndexMatrix(nEigs,dim);
for i = 1:nEigs
    c = sKernelParams.vComponentIndex(i);
    j = sKernelParams.vEigIndex(i);
    m = OneDim2MultiDimIndexMatrix(j,:);
    mPhiAnalytic(:,i) = phi(sKernelParams, c, m, mData);
    assert(~any(isnan(mPhiAnalytic(:,i))), 'Phi %d contains NaN', i);
    assert(isreal(mPhiAnalytic(:,i)), 'Phi %d is complex', i);
end


if b_normalize
    nComponents = sKernelParams.sDistParams.estNumComponents;
    mPhiAnalytic = sqrt(nComponents/nTotal)*mPhiAnalytic;
end

vLambdaAnalytic = sKernelParams.vLambdaAnalytic(1:nEigs);
fprintf('Done (took %.2f sec).\n', toc(t))
end
