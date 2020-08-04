function [mPhi_K, vLambda_K] = CalcAnalyticEigenfunctions(sParams, X)

if ~exist('X', 'var')
    x = sParams.x;
    [mX1, mX2] = meshgrid(x(:,1), x(:,2));
    X = [mX1(:) mX2(:)];    
end

if sParams.dim == 1     
    mPhi_K = zeros(length(x), sParams.PlotEigenFuncsM);
    for m = 0:sParams.PlotEigenFuncsM-1  
        mPhi_K(:,m+1) = phi_d(sParams, m, x, 1);
    end
    
    vLambda_K = zeros(sParams.PlotSpectM, 1);
    for m = 0:max(sParams.PlotSpectM,sParams.MercerM)-1
        vLambda_K(m+1) = lambda(sParams, m);
    end
    vLambda_K = vLambda_K(1:sParams.PlotSpectM);

elseif sParams.dim == 2

    mPhi_K = zeros(length(X), sParams.PlotEigenFuncsM);
    for i = 0:sParams.PlotSpectM-1
        m = OneDim2TwoDimIndex(sParams.multindexToSingleIndexMap(i+1)-1);
        mPhi_K(:,i+1) = phi(sParams, m, X);
    end
    
    vLambda_K = sParams.vLambda_K(1:sParams.PlotSpectM);
    
else
    error('cannot plot for dim > 2')
end
end
