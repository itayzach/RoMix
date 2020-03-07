function [mPhi_K, vLambda_K] = PlotAnalyticEigenfunctions(sParams, sSimParams)


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
    sgtitle('Kernel (analytic) Eigenfunctions')
    x1 = x.';
    x2 = x.';
    [mX1, mX2] = meshgrid(x1, x2);
    for m = 0:sParams.PlotEigenFuncsM-1  
        vPhi_m_x1 = phi(sParams, m, x1, 1);
        vPhi_m_x2 = phi(sParams, m, x2, 2);

        % outter product since phi(x1,x2)=phi(x1)phi(x2)
        mPhi_m_x1x2 = vPhi_m_x1.' * vPhi_m_x2; 

        subplot(2,sParams.PlotEigenFuncsM/2,m+1);
        surf(mX1, mX2, mPhi_m_x1x2, 'edgecolor', 'none')
    %         view(2)
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    end
%     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_2d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_2d.png']);
else
    error('cannot plot for dim > 2')
end
end