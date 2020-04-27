function [mPhi_K, vLambda_K] = CalcAnalyticEigenfunctions(sParams)


if sParams.dim == 1     
    mPhi_K = zeros(length(sParams.x), sParams.PlotEigenFuncsM);
    for m = 0:sParams.PlotEigenFuncsM-1  
        mPhi_K(:,m+1) = phi_d(sParams, m, sParams.x, 1);
    end
    
    vLambda_K = zeros(sParams.PlotSpectM, 1);
%     vLambda_L = zeros(sParams.PlotSpectM, 1);
    for m = 0:sParams.MercerM-1
        vLambda_K(m+1) = lambda(sParams, m);
%         vLambda_L(m+1) = lambdaL(sParams, m);
    end

elseif sParams.dim == 2
    
    vLambda_K = zeros(sParams.PlotSpectM, 1);
%     vLambda_L = zeros(sParams.PlotSpectM, 1);
    for i = 0:sParams.PlotSpectM-1
        m = OneDim2TwoDimIndex(i);
        vLambda_K(i+1) = lambda(sParams, m);
%         vLambda_L(i+1) = lambdaL(sParams, m);
    end
    [vLambda_K, sParams.multindexToSingleIndexMap] = sort(vLambda_K, 'descend');
    [mX1, mX2] = meshgrid(sParams.x(:,1), sParams.x(:,2));
    X = [mX1(:) mX2(:)];
    mPhi_K = zeros(size(sParams.x,1)*size(sParams.x,1), sParams.PlotEigenFuncsM);
    for i = 0:sParams.PlotSpectM-1
        m = OneDim2TwoDimIndex(sParams.multindexToSingleIndexMap(i+1)-1);
        mPhi_K(:,i+1) = phi(sParams, m, X);
    end

    
else
    error('cannot plot for dim > 2')
end
end
