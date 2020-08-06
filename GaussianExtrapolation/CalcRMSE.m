function vRMSE = CalcRMSE(tPhiToCompare, tPhiNumeric, compareTo)

[T, nPoints, nEigenFuncs] = size(tPhiToCompare);

mSqNorm = zeros(T, nEigenFuncs);
for t = 1:T
    for m = 1:nEigenFuncs
        if strcmp(compareTo, 'Analytic')
            mSqNorm(t,m) = norm(tPhiToCompare(t,:,m).' - sqrt(nPoints)*tPhiNumeric(t,:,m).')/sqrt(nPoints);
        elseif strcmp(compareTo, 'Nystrom')
            mSqNorm(t,m) = norm(tPhiToCompare(t,:,m).' - tPhiNumeric(t,:,m).');
        else
            error('unknown type')
        end
    end
end
vRMSE = sqrt(sum(mSqNorm,1)/T);
