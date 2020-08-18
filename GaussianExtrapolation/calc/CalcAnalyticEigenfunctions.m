function [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, mData, b_normalize)

[nTotal, dim] = size(mData);
nComponents = sKernelParams.sDistParams.estNumComponents;
mPhiAnalytic = zeros(length(mData), sSimParams.CalcEigenFuncsM);
for i = 1:sSimParams.CalcEigenFuncsM
    c = sKernelParams.vComponentIndex(i);
    j = sKernelParams.vEigIndex(i);
    if dim == 1
        m = j-1;
    elseif dim == 2
        m = OneDim2TwoDimIndex(j-1);
    end
    mPhiAnalytic(:,i) = phi(sKernelParams, c, m, mData);
end

if b_normalize
    mPhiAnalytic = sqrt(nComponents/nTotal)*mPhiAnalytic;
end

vLambdaAnalytic = sKernelParams.vLambdaAnaytic(1:sSimParams.PlotSpectM);
end
