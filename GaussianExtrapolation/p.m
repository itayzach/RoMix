function vPr = p(sParams, y, d)

if strcmp(sParams.pdf, 'gaussian')
    if exist('d', 'var')
        mu = sParams.mu(d);
        sigma = sParams.sigma(d);
    else
        mu = sParams.mu;
        sigma = sParams.sigma;
    end
    
    vPr = (1./sqrt(2*pi*sigma.^2)) .* exp( -(y-mu).^2./(2*sigma.^2) );
elseif strcmp(sParams.pdf, 'uniform')
    vPr = zeros(size(y));
    vPr(y > -0.5 & y < 0.5) = 1;
else
    error('unknown pdf')
end

end