function [mPhi_A, vLambda_A] = CalcNumericEigenvectors(sParams)
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
        elseif sParams.constsType == 3
            A = exp(-sParams.eps^2*(repmat(P',n,1) + repmat(P,1,n) - 2*(sParams.x_rand*sParams.x_rand')));
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
    
    [mPhi_A, vLambda_A] = eigs(A, sParams.PlotSpectM);
    [vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
    vLambda_A = (1/n) * vLambda_A;
    mPhi_A = sqrt(n)*mPhi_A(:,idx);

elseif sParams.dim == 2
    X = sParams.x_rand;
    P = sum(X.*X,2);
    n = length(sParams.x_rand);
    if sParams.constsType == 1
        A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(X*X'))/(2*sParams.ell^2));
    elseif sParams.constsType == 2
        A = exp(-(repmat(P',n,1) + repmat(P,1,n) - 2*(X*X'))/(2*sParams.omega^2));
    elseif sParams.constsType == 3
        A = exp(-(sParams.eps^2)*(repmat(P',n,1) + repmat(P,1,n) - 2*(X*X')));
    else
        error('Unknown constsType');
    end
    [mPhi_A, vLambda_A] = eigs(A, sParams.PlotSpectM);
    
    [vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
    vLambda_A = (1/n) * vLambda_A;
    mPhi_A = sqrt(n)*mPhi_A(:,idx);
    
else
    error('Implement for more than 2-D!');
end

end