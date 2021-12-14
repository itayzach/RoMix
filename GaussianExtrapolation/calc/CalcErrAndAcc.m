function [vRMSE, vMSE, vAcc, mErrors] = CalcErrAndAcc(tPhiToCompare, tPhiNumeric, compareTo)

[R, ~, nEigenFuncs] = size(tPhiToCompare);
mErrNorms = zeros(R, nEigenFuncs);
fprintf('*********************************************************\n');
for r = 1:R
    fprintf([compareTo '. r = %d\n'], r);
    for m = 1:nEigenFuncs
        fprintf('\t%.4f (%d)\t', max(abs(tPhiToCompare(r,:,m)))/max(abs(tPhiNumeric(r,:,m))), m-1)
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
    end
    fprintf('\n');
end
fprintf('*********************************************************\n');
vRMSE = sqrt(sum(mErrNorms.^2,1)/R).';
vMSE = (sum(mErrNorms.^2,1)/R).';
vAcc = (100*sum((1-mErrNorms),1)/R).';
mErrors = mErrNorms;
