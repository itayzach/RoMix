function arg = hermiteArg(beta,xu,mu,sigma)
arg = (1/4 + beta./2).^(1/4).*(xu-mu)./sigma;
end