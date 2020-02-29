function [mPhi_A, Lambda_A] = PlotNumericEigenvectors(sParams, sSimParams)
fig = figure;

if sParams.dim == 1
    dx = 0.1;
    x = (-5:dx:5-dx).';
    n = length(x);
%     y = x;
%     A = zeros(length(x), length(y));
%     for i = 1:length(x)
%         A(:,i) = exp(-(x(i)-y).^2./(2*sParams.l^2));
%         A(:,i) = sinc((x(i)-y)/1.5);
%         A(:,i) = (x(i)*y);
%     end

    
    % exp(-||x_i - x_j||^2 / 2\ell^2)
    P = sum(x.*x,2);
    A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(x*x'))/(2*sParams.l^2));
    
    % sinc(x_i - x_j)
%     X = repmat(x.',n,1) - repmat(x,1,n);
%     A = sinc(2*X);
    
    [mPhi_A, Lambda_A] = eig(A);
    [Lambda_A, idx] = sort(diag(Lambda_A), 'descend');
    mPhi_A = mPhi_A(:,idx);
    
    for m = 0:sSimParams.nEigenFuncsToPlot-1
        plot(x, mPhi_A(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sSimParams.nEigenFuncsToPlot - 1
            hold off
            title('Eigenvectors of (numeric) A')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
            
        end
    end
%     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_1d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_1d.png']);
%     D = diag(sum(A,2));
%     L = D - A;
%     [mPhi_L, Lambda_L] = eig(L);
%     [Lambda_L, idx] = sort(diag(Lambda_L), 'descend');
%     mPhi_L = mPhi_L(:,idx);
%     
%     fig = figure;
%     for m = 0:sSimParams.nEigenFuncsToPlot-1
%         plot(x, mPhi_L(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
%         hold on;
%         xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
%         if m == sSimParams.nEigenFuncsToPlot - 1
%             hold off
%             title('Eigenvectors of (numeric) L')
%             legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
%             print(fig, [sSimParams.outputFolder filesep 'fig1_eigenfunctions_1d'], '-depsc')
%         end
%     end
    
    
elseif sParams.dim == 2
    dx = 0.1;
    x1 = (-3:dx:3);
    x2 = (-3:dx:3);
    n = length(x1)*length(x2);
    
    [mX1, mX2] = meshgrid(x1, x2);
    X = [mX1(:) mX2(:)];
    
    P = sum(X.*X,2);
    A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(X*X'))/(2*sParams.l^2));
    
    [mPhi_A, Lambda_A] = eig(A);
    
    [Lambda_A, idx] = sort(diag(Lambda_A), 'descend');
    mPhi_A = mPhi_A(:,idx);
    
    sgtitle('Eigenvectors of (numeric) A')
    for m = 0:sSimParams.nEigenFuncsToPlot-1
        
        mPhi_Am = reshape(mPhi_A(:,m+1), length(x1), length(x2));
        
        subplot(2,sSimParams.nEigenFuncsToPlot/2,m+1);
        surf(mX1, mX2, mPhi_Am, 'edgecolor', 'none')
%         view(2)
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
%         print(fig, [sSimParams.outputFolder filesep 'fig_eigenvectors_2d'], '-depsc')
        saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_2d.png']);
    end
else
    error('Implement for more than 2-D!');
end

% fig = figure;

end