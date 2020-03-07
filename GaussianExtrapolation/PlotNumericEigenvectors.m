function [mPhi_A, vLambda_A] = PlotNumericEigenvectors(sParams, sSimParams)
fig = figure;

if sParams.dim == 1

    if strcmp(sParams.kernelType, 'sinc')
        % sinc(x_i - x_j)
        X = repmat(sParams.x_rand.',n,1) - repmat(sParams.x_rand,1,n);
        A = sinc(2*X);
    elseif strcmp(sParams.kernelType, 'exp')
        % exp(-||x_i - x_j||^2 / 2\ell^2)
        P = sum(sParams.x_rand.*sParams.x_rand,2);
        n = length(sParams.x_rand);
        if sParams.constsType == 1
            A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(sParams.x_rand*sParams.x_rand'))/(2*sParams.ell^2));
        elseif sParams.constsType == 2
            A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(sParams.x_rand*sParams.x_rand'))/(2*sParams.omega^2));
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
    
    [mPhi_A, vLambda_A] = eig(A);
    [vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
    vLambda_A = (1/length(vLambda_A)) * vLambda_A(1:sParams.PlotSpectM);
    mPhi_A = sqrt(n)*mPhi_A(:,idx);

elseif sParams.dim == 2
    %     dx = 0.1;
    %     x1 = (-3:dx:3);
    %     x2 = (-3:dx:3);
    n = 500;
    x = zeros(n, sParams.dim);
    for d = 1:sParams.dim
        if strcmp(sParams.pdf, 'gaussian')
            x(:, d) = sort((sParams.sigma(d)*randn(n, 1) + sParams.mu(d)));
        else
            error('unknown pdf')
        end
    end
    
    %     [mX1, mX2] = meshgrid(x(:, 1), x(:, 2));
    %     X = [mX1(:) mX2(:)];
    X = [x(:, 1), x(:, 2)];
    P = sum(X.*X,2);
    if sParams.constsType == 1
        A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(x*x'))/(2*sParams.ell^2));
    elseif sParams.constsType == 2
        A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(x*x'))/(2*sParams.omega^2));
    else
        error('Unknown constsType');
    end
    [mPhi_A, vLambda_A] = eig(A);
    
    [vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
    mPhi_A = mPhi_A(:,idx);
    
    sgtitle(sprintf('Eigenvectors of (numeric) A; n = %d', n))
    for m = 0:sParams.PlotEigenFuncsM-1
        
        mPhi_Am = reshape(mPhi_A(:,m+1), length(x), length(x));
        
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