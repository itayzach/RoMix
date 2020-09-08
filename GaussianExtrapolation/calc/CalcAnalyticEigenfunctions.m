function [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, mData, b_normalize)

[nTotal, dim] = size(mData);
nComponents = sKernelParams.sDistParams.estNumComponents;
mPhiAnalytic = zeros(length(mData), nEigs);
for i = 1:nEigs
    c = sKernelParams.vComponentIndex(i);
    j = sKernelParams.vEigIndex(i);
    m = OneDim2MultiDimIndex(j-1,dim);
    mPhiAnalytic(:,i) = phi(sKernelParams, c, m, mData);
end

if b_normalize
    mPhiAnalytic = sqrt(nComponents/nTotal)*mPhiAnalytic;
end

vLambdaAnalytic = sKernelParams.vLambdaAnaytic(1:nEigs);
end
