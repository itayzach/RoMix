function vPr = p(sKernelParams, y, d)

if strcmp(sKernelParams.dataDist, 'gaussian')
    if exist('d', 'var')
        mu = sKernelParams.mu(d);
        sigma = sKernelParams.sigma(d);
    else
        mu = sKernelParams.mu;
        sigma = sKernelParams.sigma;
        conv = sKernelParams.cov;
    end
    
    if sKernelParams.dim == 1
        vPr = (1./sqrt(2*pi*sigma.^2)) .* exp( -(y-mu).^2./(2*sigma.^2) );
    else
        C = 1/sqrt( ((2*pi)^sKernelParams.dim)*det(conv) );
        vPr = zeros(length(y),1);
        for i = 1:length(y)
            vPr(i) = C * exp( -0.5*(y(i,:)-mu)*(conv)^(-1)*(y(i,:)-mu).' );
        end
    end
elseif strcmp(sKernelParams.dataDist, 'uniform')
    vPr = zeros(size(y));
    vPr(y > -sKernelParams.a & y < sKernelParams.a) = 1/(2*sKernelParams.a);
%     vPr(y > 0 & y < sParams.a) = 1/sParams.a;
else
    error('unknown pdf')
end

end