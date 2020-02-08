function vPr = p(y, sigma)
vPr = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );

% vPr = prod(vPr, 2);
end