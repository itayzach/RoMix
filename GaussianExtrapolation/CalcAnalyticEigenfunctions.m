function [mPhi_K, vLambda_K] = CalcAnalyticEigenfunctions(sSimParams, sKernelParams, mData)

dim = size(mData, 2);

if dim == 1     
    mPhi_K = zeros(length(x), sSimParams.CalcEigenFuncsM);
    for m = 0:sSimParams.CalcEigenFuncsM-1  
        mPhi_K(:,m+1) = phi_d(sSimParams, m, x, 1);
    end
    
    vLambda_K = zeros(sSimParams.PlotSpectM, 1);
    for m = 0:max(sSimParams.PlotSpectM,sSimParams.MercerM)-1
        vLambda_K(m+1) = lambda(sSimParams, m);
    end
    vLambda_K = vLambda_K(1:sSimParams.PlotSpectM);

elseif dim == 2

    mPhi_K = zeros(length(mData), sSimParams.CalcEigenFuncsM);
    for i = 0:sSimParams.CalcEigenFuncsM-1
        m = OneDim2TwoDimIndex(sKernelParams.vMultindexToSingleIndexMap(i+1)-1);
        mPhi_K(:,i+1) = phi(sKernelParams, m, mData);
    end
    
    vLambda_K = sKernelParams.vLambdaAnaytic(1:sSimParams.PlotSpectM);
    
else
    error('cannot plot for dim > 2')
end


end
