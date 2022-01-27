function [lambda_m] = lambdaD(sKernelParams, c, m)

beta = sKernelParams.beta{c};
vLambda_m = sqrt(2./(1+beta+sqrt(1+2*beta))) .* (beta./(1+beta+sqrt(1+2*beta))).^m;
lambda_m = prod(vLambda_m);

end