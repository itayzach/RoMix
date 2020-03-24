function vPr = p(sParams, y, d)

if strcmp(sParams.dataDist, 'gaussian')
    if exist('d', 'var')
        mu = sParams.mu(d);
        sigma = sParams.sigma(d);
    else
        mu = sParams.mu;
        sigma = sParams.sigma;
    end
    
%     if sParams.dim == 1
        vPr = (1./sqrt(2*pi*sigma.^2)) .* exp( -(y-mu).^2./(2*sigma.^2) );
%     else
%         C = 1/sqrt( ((2*pi)^sParams.dim)*det(sigma) );
%         vPr = C * exp( -0.5*(y-mu)*((sigma)\(y-mu).') );
%     end
elseif strcmp(sParams.dataDist, 'uniform')
    vPr = zeros(size(y));
    vPr(y > -0.5 & y < 0.5) = 1;
else
    error('unknown pdf')
end

end