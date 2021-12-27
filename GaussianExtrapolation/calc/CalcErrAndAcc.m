function [vRMSE, vMSE, vAcc, mErrors, vCoh] = CalcErrAndAcc(tPhiToCompare, tPhiNumeric, compareTo)
if ismatrix(tPhiToCompare)
    tPhiToCompare = reshape(tPhiToCompare,1,size(tPhiToCompare,1),size(tPhiToCompare,2));
    tPhiNumeric = reshape(tPhiNumeric,1,size(tPhiNumeric,1),size(tPhiNumeric,2));
end
[R, ~, nEigenFuncs] = size(tPhiToCompare);
mErrNormed = zeros(R, nEigenFuncs);
mCoherence = zeros(R, nEigenFuncs);
fprintf('*********************************************************\n');
for r = 1:R
    fprintf(['Accuracy of ', compareTo '. r = %d\n'], r);
    mPhiToCompare = squeeze(tPhiToCompare(r,:,:));
    mPhiNumeric = squeeze(tPhiNumeric(r,:,:));
    mErrNormed(r,:) = vecnorm(mPhiToCompare-mPhiNumeric)./vecnorm(mPhiToCompare);
    mCoherence(r,:) = (vecnorm(mPhiToCompare-mPhiNumeric).^2)./(vecnorm(mPhiToCompare).*vecnorm(mPhiNumeric));
    for m = 1:nEigenFuncs
        fprintf('\t(m=%2d) %8.4f%%\t', m-1, 100*(1-mErrNormed(r,m)))
        if (mod(m,8) == 0) && m < nEigenFuncs
            fprintf('\n');
        end
    end
    fprintf('\n');
end
fprintf('*********************************************************\n');

% The sum is over the iterations, r=1:R
vRMSE = sqrt(sum(mErrNormed.^2,1)/R).';
vMSE = (sum(mErrNormed.^2,1)/R).';
vCoh = (sum(mCoherence,1)/R).';
vAcc = (100*sum((1-mErrNormed),1)/R).';
mErrors = mErrNormed.';
end


