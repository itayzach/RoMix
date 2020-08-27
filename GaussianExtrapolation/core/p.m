function vDensity = p(sDistParams, c, y, d)

if strcmp(sDistParams.estDataDist, 'Gaussian')
    if exist('d', 'var')
        mu = sDistParams.mu{c}(d);
        sigma = sDistParams.sigma{c}(d);
    else
        mu = sDistParams.mu{c};
        sigma = sDistParams.sigma{c};
        cov = sDistParams.cov{c};
    end
    
    if sDistParams.dim == 1
        vDensity = (1./sqrt(2*pi*sigma.^2)) .* exp( -(y-mu).^2./(2*sigma.^2) );
    else
        C = 1/sqrt( ((2*pi)^sDistParams.dim)*det(cov) );
        vDensity = zeros(length(y),1);
        for i = 1:length(y)
            vDensity(i) = C * exp( -0.5*(y(i,:)-mu)*(cov)^(-1)*(y(i,:)-mu).' );
        end
    end
elseif strcmp(sDistParams.dataDist, 'uniform')
    vDensity = zeros(size(y));
    vDensity(y > -sDistParams.a & y < sDistParams.a) = 1/(2*sDistParams.a);
%     vPr(y > 0 & y < sParams.a) = 1/sParams.a;
else
    error('unknown pdf')
end

end