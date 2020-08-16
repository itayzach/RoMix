function [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, mData)

[nTotal, dim] = size(mData);
nComponents = sKernelParams.sDistParams.estNumComponents;
mPhiAnalytic = zeros(length(mData), sSimParams.CalcEigenFuncsM);
for i = 1:nComponents*sSimParams.CalcEigenFuncsM
    if dim == 1
        m = i-1;
    elseif dim == 2
        m = OneDim2TwoDimIndex(sKernelParams.vMultindexToSingleIndexMap(i)-1);
    end
    mPhiAnalytic(:,i) = sqrt(nComponents/nTotal)*phi(sKernelParams, m, mData);
end   

vLambdaAnalytic = sKernelParams.vLambdaAnaytic(1:sSimParams.PlotSpectM);
end
