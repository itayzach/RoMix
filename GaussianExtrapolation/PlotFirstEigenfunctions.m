function [] = PlotFirstEigenfunctions(sParams, sSimParams)
dx = 0.01;
x = (-5:dx:5-dx).';
fig1 = figure(1);

pltIdx = 1;
for m = 0:sSimParams.nEigenFuncsToPlot-1  
    if sParams.dim == 1  
        if m == 0
            vP_x = p(x, sParams.sigma);
            
            hold on
            % subplot(floor(sSimParams.nEigenFuncsToPlot/2), floor(sSimParams.nEigenFuncsToPlot/2)+1, pltIdx)
            plot(x, vP_x, '-.', 'LineWidth', 2, 'DisplayName', '$p(x)$');
            pltIdx = pltIdx + 1;
        end   
        
        [vPhi_m_x, ~] = phi(sParams.a, sParams.b, m, x);
        % subplot(floor(sSimParams.nEigenFuncsToPlot/2), floor(sSimParams.nEigenFuncsToPlot/2)+1, pltIdx)
        plot(x, vPhi_m_x, 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sSimParams.nEigenFuncsToPlot - 1
            hold off
            %     title('Eigenfunctions')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
            print(fig1, [sSimParams.outputFolder filesep 'fig1_eigenfunctions'], '-depsc')
        end
        pltIdx = pltIdx + 1;
    elseif sParams.dim == 2
        x1 = x.';
        x2 = x.';
        [mX1, mX2] = meshgrid(x1, x2);
        [vPhi_m_x1, ~] = phi(sParams.a, sParams.b, m, x1);
        [vPhi_m_x2, ~] = phi(sParams.a, sParams.b, m, x2);
        
        % outter product since phi(x1,x2)=phi(x1)phi(x2)
        mPhi_m_x1x2 = vPhi_m_x1.' * vPhi_m_x2; 
        
        subplot(2,2,m+1);
        surf(mX1, mX2, mPhi_m_x1x2, 'edgecolor', 'none')
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    else
        error('cannot plot for dim > 2')
    end
    
end
end