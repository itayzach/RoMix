function [mPhiAnalytic, vLambdaAnalytic, t] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, mData)

[nTotal, dim] = size(mData);
OneDim2MultiDimIndexMatrix = LoadOneDim2MultiDimIndexMatrix(nEigs,dim);
fprintf('Calculating %d eigenfunctions (max m = %d) d = %d, n = %d... ',nEigs, max(OneDim2MultiDimIndexMatrix(:)), dim, nTotal)
ts = tic;
mPhiAnalytic = zeros(nTotal, nEigs);
for i = 1:nEigs
    c = sKernelParams.vComponentIndex(i);
    j = sKernelParams.vEigIndex(i);
    m = OneDim2MultiDimIndexMatrix(j,:);
    mPhiAnalytic(:,i) = phiD(sKernelParams, c, m, mData);
%     isalmostequal(mPhiAnalytic(:,i),phi(sKernelParams, c, m, mData),1e-13)
    assert(~any(isnan(mPhiAnalytic(:,i))), 'Phi %d contains NaN', i);
    assert(isreal(mPhiAnalytic(:,i)), 'Phi %d is complex', i);
end

nComponents = sKernelParams.sDistParams.estNumComponents;
mPhiAnalytic = sqrt(nComponents)*mPhiAnalytic;

vLambdaAnalytic = sKernelParams.vLambdaAnalytic(1:nEigs);
t = toc(ts);
fprintf('Done (took %.2f sec).\n', t)
end
