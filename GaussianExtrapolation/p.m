function vPr = p(sParams, y, d)

if strcmp(sParams.dataDist, 'gaussian')
    if exist('d', 'var')
        mu = sParams.mu(d);
        sigma = sParams.sigma(d);
    else
        mu = sParams.mu;
        sigma = sParams.sigma;
        conv = sParams.cov;
    end
    
    if sParams.dim == 1
        vPr = (1./sqrt(2*pi*sigma.^2)) .* exp( -(y-mu).^2./(2*sigma.^2) );
    else
        C = 1/sqrt( ((2*pi)^sParams.dim)*det(conv) );
        vPr = zeros(length(y),1);
        for i = 1:length(y)
            vPr(i) = C * exp( -0.5*(y(i,:)-mu)*(conv)^(-1)*(y(i,:)-mu).' );
        end
    end
elseif strcmp(sParams.dataDist, 'uniform')
    vPr = zeros(size(y));
    vPr(y > -sParams.a & y < sParams.a) = 1/(2*sParams.a);
%     vPr(y > 0 & y < sParams.a) = 1/sParams.a;
else
    error('unknown pdf')
end

end