function [mPhi_K, vLambda_K] = CalcAnalyticEigenfunctions(sParams)


if sParams.dim == 1     
    mPhi_K = zeros(length(sParams.x), sParams.PlotEigenFuncsM);
    for m = 0:sParams.PlotEigenFuncsM-1  
        mPhi_K(:,m+1) = phi_d(sParams, m, sParams.x, 1);
    end
    
    vLambda_K = zeros(sParams.PlotSpectM, 1);
    for m = 0:sParams.PlotSpectM-1
        vLambda_K(m+1) = prod(lambda(sParams, m), 2);
    end

elseif sParams.dim == 2
    
    vLambda_K = zeros(sParams.PlotSpectM, 1);
    for i = 0:sParams.PlotSpectM-1
        m = OneDim2TwoDimIndex(i, sParams.dim);
        vLambda_K(i+1) = lambda(sParams, m);
    end
    
    [mX1, mX2] = meshgrid(sParams.x(:,1), sParams.x(:,2));
    X = [mX1(:) mX2(:)];
    mPhi_K = zeros(size(sParams.x,1)*size(sParams.x,1), sParams.PlotEigenFuncsM);
    for i = 0:sParams.PlotSpectM-1
        m = OneDim2TwoDimIndex(i, sParams.dim);
        mPhi_K(:,i+1) = phi(sParams, m, X);
    end

    
else
    error('cannot plot for dim > 2')
end
end