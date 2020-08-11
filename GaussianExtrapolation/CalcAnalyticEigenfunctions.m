function [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, mData)

[nTotal, dim] = size(mData);

mPhiAnalytic = zeros(length(mData), sSimParams.CalcEigenFuncsM);
for i = 0:sSimParams.CalcEigenFuncsM-1
    if dim == 1
        m = i;
    elseif dim == 2
        m = OneDim2TwoDimIndex(sKernelParams.vMultindexToSingleIndexMap(i+1)-1);
    end
    mPhiAnalytic(:,i+1) = (1/sqrt(nTotal))*phi(sKernelParams, m, mData);
end   

vLambdaAnalytic = sKernelParams.vLambdaAnaytic(1:sSimParams.PlotSpectM);
end
