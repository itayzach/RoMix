function x = GenerateGaussianData(dim, nComponents, n)
if dim >= 3
    sigma = 1;
    mu = 0;
    x = sigma*randn(n,dim) + mu;
    
elseif dim == 2
    if nComponents == 1
        cov = [0.25    0.01;
            0.01   0.25];
        mu  = [0; 0];
        x = mvnrnd(mu, cov, n);
    elseif nComponents == 2
        cov1 = [0.25    0.01;
            0.01   0.5];
        mu1  = [0; 0];
        x1 = mvnrnd(mu1, cov1, n);
        
        cov2 = [0.7    -0.1;
            -0.1   0.1];
        mu2  = [5; 5];
        x2 = mvnrnd(mu2, cov2, n);
        
        vSel = rand(n, 1) < 0.5;
        x = (1-vSel).*x1 + vSel.*x2;
        
    else
        error('not supported')
    end
elseif dim == 1
    if nComponents == 1
        sigma = 1;
        mu = 0;
        x = sigma*randn(n, 1) + mu;
    elseif nComponents == 2
        mu = [2; 7];       % Means
        sigma = [0.4 0.5]; % Covariances
        vSel = rand(n, 1) < 0.5;
        x = (1-vSel).*(sigma(1)*randn(n, 1) + mu(1)) + ...
            vSel.*(sigma(2)*randn(n, 1) + mu(2));
    else
        error('not supported')
    end
else
    error('generate random data for more than 2D')
end
end