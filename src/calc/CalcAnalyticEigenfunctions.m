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
end

vNormFactor = reshape(sKernelParams.sDistParams.GMModel.ComponentProportion(sKernelParams.vComponentIndex), 1, []);
mPhiAnalytic = 1./sqrt(vNormFactor).*mPhiAnalytic;
vLambdaAnalytic = sKernelParams.vLambdaAnalytic;
t = toc(ts);

assert(~any(isnan(mPhiAnalytic(:))), 'Phi contains NaN');
assert(isreal(mPhiAnalytic(:)), 'Phi is complex');
fprintf('Done (took %.2f sec).\n', t)
end