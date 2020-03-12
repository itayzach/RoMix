function [mPhi_K, vLambda_K] = CalcAnalyticEigenfunctions(sParams, sSimParams)


mPhi_K = zeros(length(sParams.x), sParams.PlotEigenFuncsM);
vLambda_K = zeros(sParams.PlotEigenFuncsM, 1);

if sParams.dim == 1     
    for m = 0:sParams.PlotEigenFuncsM-1  
        mPhi_K(:,m+1) = phi(sParams, m, sParams.x, 1);
    end
    vLambda_K = zeros(sParams.PlotSpectM, 1);
    for m = 0:sParams.PlotSpectM-1
        vLambda_K(m+1) = prod(lambda(sParams, m), 2);
    end

elseif sParams.dim == 2
    
else
    error('cannot plot for dim > 2')
end
end