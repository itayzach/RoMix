function [vRMSE, vMSE, vAcc, mErrors, vCoh, vAccStd] = CalcErrAndAcc(tPhiToCompare, tPhiNumeric, compareTo)
[n, nEigenFuncs, R] = size(tPhiToCompare);
if isequal(tPhiToCompare(:,1,1),floor(tPhiToCompare(:,1,1)))
    errFunc = '0-1';
else
    errFunc = 'norm';
end

mErrNormed = zeros(R, nEigenFuncs);
mCoherence = zeros(R, nEigenFuncs);
fprintf('*********************************************************\n');
for r = 1:R
    fprintf(['Accuracy of ', compareTo '. r = %d\n'], r);
    mPhiToCompare = squeeze(tPhiToCompare(:,:,r));
    mPhiNumeric = squeeze(tPhiNumeric(:,:,r));
    if strcmp(errFunc, 'norm')
        mErrNormed(r,:) = vecnorm(mPhiToCompare-mPhiNumeric,2)./vecnorm(mPhiToCompare,2);
        mCoherence(r,:) = (vecnorm(mPhiToCompare-mPhiNumeric,2).^2)./(vecnorm(mPhiToCompare,2).*vecnorm(mPhiNumeric,2));
    elseif strcmp(errFunc, '0-1')
        mErrNormed(r,:) = sum(mPhiToCompare ~= mPhiNumeric) / n;
        mCoherence(r,:) = sum(mPhiToCompare ~= mPhiNumeric) / n;
    end
    for m = 1:nEigenFuncs
        fprintf('\t(m=%2d) %6.4f%%\t', m-1, 100*(1-mErrNormed(r,m)))
        if (mod(m,8) == 0) && m < nEigenFuncs
            fprintf('\n');
        end
    end
    fprintf('\n');
end
% fprintf('*********************************************************\n');

% The sum is over the iterations, r=1:R
vRMSE = sqrt(sum(mErrNormed.^2,1)/R).';
vMSE = (sum(mErrNormed.^2,1)/R).';
vCoh = (sum(mCoherence,1)/R).';
vAcc = (100*sum((1-mErrNormed),1)/R).';
mErrors = mErrNormed.';
vAccStd = 100*std(mErrNormed,[],1).';
end


