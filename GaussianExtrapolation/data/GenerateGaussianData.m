function x = GenerateGaussianData(dim, nComponents, n, mu, sigma)
if dim >= 3
    assert(nComponents == 1, 'not implemented more than 1 component')
    x = mvnrnd(mu, sigma, n);
elseif dim == 2
    if nComponents == 1
        x = mvnrnd(mu{1}, sigma{1}, n);
    elseif nComponents == 2
        msgBoxMsg = 'Note that cov and mu are fixed in this function';
        msgBoxTitle = 'GenerateGaussianData';
        msgBoxIcon = 'warn';
        uiwait(msgbox(msgBoxMsg, msgBoxTitle, msgBoxIcon))
        cov1 = [0.25,    0.01;
                0.01,    0.5];
        mu1  = [0; 0];
        x1 = mvnrnd(mu1, cov1, n);
        
        cov2 = [0.7,    -0.1;
               -0.1,     0.1];
        mu2  = [5; 5];
        x2 = mvnrnd(mu2, cov2, n);
        
        vSel = rand(n, 1) < 0.5;
        x = (1-vSel).*x1 + vSel.*x2;
        
    else
        error('not supported')
    end
elseif dim == 1
    x = zeros(n,1);
    nPtsPerComp = round(n/nComponents);
    for c = 1:nComponents
        vPtsInd = (1:nPtsPerComp) + (c-1)*nPtsPerComp;
        x(vPtsInd) = sigma{c}*randn(nPtsPerComp, 1) + mu{c};
    end
    vRandInd = randperm(n);
    x = x(vRandInd);
else
    error('generate random data for more than 3D')
end
end