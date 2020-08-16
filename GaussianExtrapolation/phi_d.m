function [vPhi_m] = phi_d(sKernelParams, m, x, d)
if ~isvector(x)
    error('x has to be a vector');
end

if strcmp(sKernelParams.kernelType, 'gaussian')
    if sKernelParams.constsType == 1
        a = sKernelParams.a(d);
        b = sKernelParams.b;

        % Calculate parameters
        c = sqrt(a^2 + 2*a*b);

        % m-th eigenfunction
        if m <= 170
            % factorial(171) = Inf, but hermite() is more efficient than hermiteH...
            vHm = hermite(m, sqrt(2*c)*x);
        else
            warning('m > 170, using hermiteH instead of factorial hermite')
            vHm = hermiteH(m, sqrt(2*c)*x);
        end

        normFactor = (1/sqrt(2^m*factorial(m)*sqrt(a/c)));
        vPhi_m = normFactor * exp( -(c-a)*x.^2 ) .* vHm;

    elseif sKernelParams.constsType == 2
        if exist('d', 'var')
            mu = sKernelParams.sDistParams.mu_1D(d);
            sigma = sKernelParams.sDistParams.sigma(d);
            beta = sKernelParams.beta(d);
        else
            mu = sKernelParams.sDistParams.mu_1D;
            sigma = sKernelParams.sDistParams.sigma;
            beta = sKernelParams.beta;
        end

        normFactor = (1+2*beta)^(1/8)/sqrt(2^m*factorial(m));
        vHm = hermite(m, (1/4 + beta/2)^(1/4)*(x-mu)/sigma);
        vPhi_m = normFactor * exp( -((x-mu).^2/(2*sigma^2)) * ((sqrt(1+2*beta)-1)/2) ) .* vHm;
    elseif sKernelParams.constsType == 3
        if exist('d', 'var')
            mu = sKernelParams.mu(d);
            sigma = sKernelParams.sDistParams.sigma(d);     
            alpha = sKernelParams.alpha(d); % Should be here, or alpha(d)?
        else
            mu = sKernelParams.mu;
            sigma = sKernelParams.sDistParams.sigma;
        end    
        eps = sKernelParams.eps;


        normFactor = (1+(2*eps/alpha)^2)^(1/8) / sqrt(2^m*factorial(m));
        vHm = hermite(m, (1+(2*eps/alpha)^2).^(1/4).*alpha.*x );
        vPhi_m = normFactor * exp( -(sqrt(1+(2*eps/alpha)^2)-1)*alpha^2*x.^2/2 ) .* vHm;
    else
        error('Unknown constsType');
    end
elseif strcmp(sKernelParams.kernelType, 'sinc')
%     vPhi_m = real(sqrt(pi./(2*sParams.a*x)) .* besselj(m+0.5, sParams.a*x));
%     vPhi_m(x < 0) = -vPhi_m(x < 0);
%     if m == 0
%         vPhi_m(x == 0) = 1;
%     else
%         vPhi_m(x == 0) = 0;
%     end
    vPhi_m = besselj(m, sKernelParams.a*x);
else
    error('unknown kernelType');
end
end