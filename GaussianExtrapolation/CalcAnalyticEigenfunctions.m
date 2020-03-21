function [mPhi_K, vLambda_K] = CalcAnalyticEigenfunctions(sParams, sSimParams)


if sParams.dim == 1     
    mPhi_K = zeros(length(sParams.x), sParams.PlotEigenFuncsM);
    for m = 0:sParams.PlotEigenFuncsM-1  
        mPhi_K(:,m+1) = phi(sParams, m, sParams.x, 1);
    end
    
    vLambda_K = zeros(sParams.PlotSpectM, 1);
    for m = 0:sParams.PlotSpectM-1
        vLambda_K(m+1) = prod(lambda(sParams, m), 2);
    end

elseif sParams.dim == 2
    
    vLambda_Kd = zeros(sParams.PlotSpectM, sParams.dim);
    for m = 0:sParams.PlotSpectM-1
        vLambda_Kd(m+1,:) = lambda(sParams, m);
    end
    
    vLambda_K = vLambda_Kd(:,1)*vLambda_Kd(:,1).';
    [vLambda_K, idx] = sort(vLambda_K(:), 'descend');
    vLambda_K = vLambda_K(1:sParams.PlotSpectM);
    
    tPhi_x = zeros(size(sParams.x,1), sParams.PlotEigenFuncsM, sParams.dim);
    for m = 0:sParams.PlotEigenFuncsM-1 
        for d = 1:sParams.dim
            tPhi_x(:,m+1,d) = phi(sParams, m, sParams.x(:,d), d);
            assert(~any(isnan(squeeze(tPhi_x(:,m+1,d)))));
        end
    end
    
    mPhi_K = zeros(size(sParams.x,1)*size(sParams.x,1), sParams.PlotEigenFuncsM);  
    mPhi_m_x = repmat(tPhi_x(:,m+1,1), 1, size(sParams.x,1));
    mPhi_m_y = repmat(tPhi_x(:,m+1,2).', size(sParams.x,1), 1);
    for m = 0:sParams.PlotEigenFuncsM-1 
        
        mPhi_m_xy = mPhi_m_x.*mPhi_m_y;
        mPhi_K(:,m+1) = mPhi_m_xy(:);
        
%         for i = 1:size(sParams.x,1)
%             for j = 1:size(sParams.x,1)
%                 mPhi_K2_tmp(i,j) = tPhi_x(i,m+1,1)*tPhi_x(j,m+1,2);
%             end
%         end
%         mPhi_K2(:,m+1) = mPhi_K2_tmp(:);
        
    end
%     mPhi_K(:,m+1) = prod(tPhi_xx(:,m+1,:),3);

    
else
    error('cannot plot for dim > 2')
end
end