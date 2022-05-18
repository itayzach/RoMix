function x = GenerateGaussianData(dim, nComponents, n, mu, sigma, p)
mMu = cell2mat(mu');
mCov = reshape(cell2mat(sigma), dim, dim, nComponents);
mP = cell2mat(p);
gm = gmdistribution(mMu, mCov, mP);
x = random(gm, n);
end