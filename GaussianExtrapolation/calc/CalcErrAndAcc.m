function [vRMSE, vMSE, vAcc, mErrors] = CalcErrAndAcc(tPhiToCompare, tPhiNumeric, compareTo)
if ismatrix(tPhiToCompare)
    tPhiToCompare = reshape(tPhiToCompare,1,size(tPhiToCompare,1),size(tPhiToCompare,2));
    tPhiNumeric = reshape(tPhiNumeric,1,size(tPhiNumeric,1),size(tPhiNumeric,2));
end
[R, ~, nEigenFuncs] = size(tPhiToCompare);
mErrNorms = zeros(R, nEigenFuncs);
fprintf('*********************************************************\n');
for r = 1:R
    fprintf(['Accuracy of ', compareTo '. r = %d\n'], r);
    for m = 1:nEigenFuncs
%         fprintf('\t%.4f (%d)\t', max(abs(tPhiToCompare(r,:,m)))/max(abs(tPhiNumeric(r,:,m))), m-1)
        if (mod(m,10) == 0) && m < nEigenFuncs
            fprintf('\n');
        end
%         if strcmp(compareTo, 'Analytic')
%             mSqNorm(r,m) = sqrt((1/nPoints)*sum((tPhiToCompare(r,:,m).' - tPhiNumeric(r,:,m).').^2));
%         elseif strcmp(compareTo, 'Nystrom')
%             mSqNorm(r,m) = sqrt((1/nPoints)*sum((tPhiToCompare(r,:,m).' - tPhiNumeric(r,:,m).').^2));
%         else
%             error('unknown type')
%         end
        mErrNorms(r,m) = min(norm(tPhiToCompare(r,:,m) - tPhiNumeric(r,:,m)),...
                           norm(tPhiToCompare(r,:,m) + tPhiNumeric(r,:,m)))/norm(tPhiToCompare(r,:,m));
        fprintf('\t(m=%d) %.4f%%\t', m-1, 100*(1-mErrNorms(r,m)))
    end
    fprintf('\n');
end
fprintf('*********************************************************\n');
vRMSE = sqrt(sum(mErrNorms.^2,1)/R).';
vMSE = (sum(mErrNorms.^2,1)/R).';
vAcc = (100*sum((1-mErrNorms),1)/R).';
mErrors = mErrNorms;
