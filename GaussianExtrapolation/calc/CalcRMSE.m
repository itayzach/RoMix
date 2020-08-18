function vRMSE = CalcRMSE(tPhiToCompare, tPhiNumeric, compareTo)

[T, nPoints, nEigenFuncs] = size(tPhiToCompare);

mSqNorm = zeros(T, nEigenFuncs);
fprintf('*********************************************************\n');
for t = 1:T
    fprintf([compareTo '. t = %d\n'], t);
    for m = 1:nEigenFuncs
        fprintf('m = %d; Ratio: %f\n', m-1, max(tPhiToCompare(t,:,m))/max(tPhiNumeric(t,:,m)))
        if strcmp(compareTo, 'Analytic')
            mSqNorm(t,m) = sqrt((1/nPoints)*sum((tPhiToCompare(t,:,m).' - tPhiNumeric(t,:,m).').^2));
        elseif strcmp(compareTo, 'Nystrom')
            mSqNorm(t,m) = sqrt((1/nPoints)*sum((tPhiToCompare(t,:,m).' - tPhiNumeric(t,:,m).').^2));
        else
            error('unknown type')
        end
    end
    fprintf('\n');
end
fprintf('*********************************************************\n');
vRMSE = sqrt(sum(mSqNorm,1)/T);
