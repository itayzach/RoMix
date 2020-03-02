function vPr = p(sParams, y)
mu = sParams.mu;
sigma = sParams.sigma;

vPr = (1/sqrt(2*pi*sigma^2)) * exp( -(y-mu).^2/(2*sigma^2) );
end