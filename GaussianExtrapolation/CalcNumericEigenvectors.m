function [mPhi_A, vLambda_A] = CalcNumericEigenvectors(sParams)
x = sParams.x_rand;
n = length(x);

if strcmp(sParams.kernelType, 'sinc')
    % sinc(x_i - x_j)
    dist = pdist2(x, x);
%     X = repmat(x.',n,1) - repmat(x,1,n);
    A = sin(sParams.a*dist)./(pi*dist);
    A(dist == 0) = sParams.a;
elseif strcmp(sParams.kernelType, 'exp')
    % exp(-||x_i - x_j||^2 / 2\ell^2)
    dist = pdist2(x, x);
    if sParams.constsType == 1
        A = exp(-dist.^2/(2*sParams.ell^2));
    elseif sParams.constsType == 2
        A = exp(-dist.^2/(2*sParams.omega^2));
    elseif sParams.constsType == 3
        A = exp(-sParams.eps^2 * dist.^2);
    else
        error('Unknown constsType');
    end
% elseif strcmp(sParams.kernelType, 'brownian')
%     x = sParams.x_rand;
%     n = length(x);
%     dist = pdist2(x, x);
%     A = sParams.a^2 - x*x.' -0.5*dist;
else
    error('Unknown kernelType')
end

[mPhi_A, vLambda_A] = eigs(A, sParams.PlotSpectM);
[vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
vLambda_A = (1/n) * vLambda_A;
mPhi_A = sqrt(n)*mPhi_A(:,idx);

end