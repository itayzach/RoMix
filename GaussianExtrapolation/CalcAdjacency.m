function A = CalcAdjacency(sParams, x, y)

if ~exist('y', 'var')
    y = x;
end

if strcmp(sParams.kernelType, 'sinc')
    % sinc(x_i - x_j)
    dist = pdist2(x, y);
    A = sin(sParams.a*dist)./(pi*dist);
    A(dist == 0) = sParams.a;
elseif strcmp(sParams.kernelType, 'exp')
    % exp(-||x_i - x_j||^2 / 2\ell^2)
    dist = pdist2(x, y);
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
end