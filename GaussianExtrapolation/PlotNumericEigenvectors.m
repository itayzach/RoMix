function [mPhi_A, Lambda_A_avg] = PlotNumericEigenvectors(sParams, sSimParams)
fig = figure;

if sParams.dim == 1
    nExp = 1;
    n = 5000;
    Lambda_Ai = zeros(n, nExp);
    for i = 1:nExp
        x = zeros(n, sParams.dim);
        for d = 1:sParams.dim
            if strcmp(sParams.pdf, 'gaussian')
                x(:, d) = sort((sParams.sigma*randn(n, 1) + sParams.mu));
                xMin = min(x);
                xMax = max(x);
            elseif strcmp(sParams.pdf, 'uniform')
                xMin = -0.5;
                xMax = 0.5;
                x(:, d) = (xMax - xMin)*sort(rand(n, 1)) + xMin;
            else
                error('unknown pdf')
            end
        end
        
        if strcmp(sParams.kernelType, 'sinc')
            % sinc(x_i - x_j)
            X = repmat(x.',n,1) - repmat(x,1,n);
            A = sinc(2*X);
        elseif strcmp(sParams.kernelType, 'exp')
            % exp(-||x_i - x_j||^2 / 2\ell^2)
            P = sum(x.*x,2);
            if sParams.constsType == 1
                A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(x*x'))/(2*sParams.l^2));
            elseif sParams.constsType == 2
                A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(x*x'))/(2*sParams.omega^2));
            else
                error('Unknown constsType');
            end
        else
            error('Unknown kernelType')
        end
        % Normalize
    %     D1 = sum(A, 1);
    %     D2 = sum(A, 2);
    %     D = sqrt(D2*D1);
    %     A_normalized = A./D;
    %     A_normalized = diag(D1)^(-1/2) * A * diag(D1)^(-1/2);
    %     [mPhi_A, Lambda_A] = eig(A_normalized);

        [mPhi_A, Lambda_A] = eig(A);
        [Lambda_A, idx] = sort(diag(Lambda_A), 'descend');
        Lambda_Ai(:,i) = Lambda_A;
        mPhi_A = sqrt(n)*mPhi_A(:,idx);
        
        for m = 0:sParams.PlotEigenFuncsM-1            
            plot(x, mPhi_A(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
            hold on
            xlim([xMin xMax]);
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
            if m == sParams.PlotEigenFuncsM - 1
                histogram(x, 'Normalization', 'pdf', 'LineStyle', ':', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', '$\hat{p}(x)$');
                hold off
                title('Eigenvectors of (numeric) A')
                legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
            end
        end
    end
    Lambda_A_avg = mean(Lambda_Ai, 2);
%     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_1d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_1d.png']);
    
    
elseif sParams.dim == 2
    dx = 0.1;
    x1 = (-3:dx:3);
    x2 = (-3:dx:3);
    n = length(x1)*length(x2);
    
    [mX1, mX2] = meshgrid(x1, x2);
    X = [mX1(:) mX2(:)];
    
    P = sum(X.*X,2);
    A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(X*X'))/(2*sParams.omega^2));
    
    [mPhi_A, Lambda_A] = eig(A);
    
    [Lambda_A, idx] = sort(diag(Lambda_A), 'descend');
    mPhi_A = mPhi_A(:,idx);
    
    sgtitle(sprintf('Eigenvectors of (numeric) A; n = %d', n))
    for m = 0:sParams.PlotEigenFuncsM-1
        
        mPhi_Am = reshape(mPhi_A(:,m+1), length(x1), length(x2));
        
        subplot(2,sParams.PlotEigenFuncsM/2,m+1);
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