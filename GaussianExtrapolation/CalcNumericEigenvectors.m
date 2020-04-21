function [mPhi_A, vLambda_A] = CalcNumericEigenvectors(sParams)

n = length(sParams.x_rand);
if strcmp(sParams.kernelType, 'sinc')
    % sinc(x_i - x_j)
    X = repmat(sParams.x_rand.',n,1) - repmat(sParams.x_rand,1,n);
    A = sinc(2*X);
elseif strcmp(sParams.kernelType, 'exp')
    % exp(-||x_i - x_j||^2 / 2\ell^2)
    dist = pdist2(sParams.x_rand, sParams.x_rand);
    if sParams.constsType == 1
        A = exp(-dist.^2/(2*sParams.ell^2));
    elseif sParams.constsType == 2
        A = exp(-dist.^2/(2*sParams.omega^2));
    elseif sParams.constsType == 3
        A = exp(-sParams.eps^2 * dist.^2);
    else
        error('Unknown constsType');
    end
else
    error('Unknown kernelType')
end

[mPhi_A, vLambda_A] = eigs(A, sParams.PlotSpectM);
[vLambda_A, idx] = sort(diag(vLambda_A), 'descend');
vLambda_A = (1/n) * vLambda_A;
mPhi_A = sqrt(n)*mPhi_A(:,idx);

end