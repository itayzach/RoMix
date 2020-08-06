function [mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, mData)

dim = size(mData, 2);

if dim == 1     
    mPhiAnalytic = zeros(length(x), sSimParams.CalcEigenFuncsM);
    for m = 0:sSimParams.CalcEigenFuncsM-1  
        mPhiAnalytic(:,m+1) = phi_d(sSimParams, m, x, 1);
    end
    
    vLambdaAnalytic = zeros(sSimParams.PlotSpectM, 1);
    for m = 0:max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM)-1
        vLambdaAnalytic(m+1) = lambda(sSimParams, m);
    end
    vLambdaAnalytic = vLambdaAnalytic(1:sSimParams.PlotSpectM);

elseif dim == 2

    mPhiAnalytic = zeros(length(mData), sSimParams.CalcEigenFuncsM);
    for i = 0:sSimParams.CalcEigenFuncsM-1
        m = OneDim2TwoDimIndex(sKernelParams.vMultindexToSingleIndexMap(i+1)-1);
        mPhiAnalytic(:,i+1) = phi(sKernelParams, m, mData);
    end
    
    vLambdaAnalytic = sKernelParams.vLambdaAnaytic(1:sSimParams.PlotSpectM);
    
else
    error('cannot plot for dim > 2')
end


end
