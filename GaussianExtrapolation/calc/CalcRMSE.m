function vRMSE = CalcRMSE(tPhiToCompare, tPhiNumeric, compareTo)

% [R, ~, nEigenFuncs] = size(tPhiToCompare);
% 
% mSqNorm = zeros(R, nEigenFuncs);
% fprintf('*********************************************************\n');
% for r = 1:R
%     fprintf([compareTo '. r = %d\n'], r);
%     for m = 1:nEigenFuncs
%         fprintf('m = %d; Ratio: %f\n', m-1, max(tPhiToCompare(r,:,m))/max(tPhiNumeric(r,:,m)))
%         mSqNorm(r,m) = norm(tPhiToCompare(r,:,m) - tPhiNumeric(r,:,m));
%     end
%     fprintf('\n');
% end
% fprintf('*********************************************************\n');
% vRMSE = sqrt(sum(mSqNorm.^2,1)/R);


[R, nPoints, nEigenFuncs] = size(tPhiToCompare);

mSqNorm = zeros(R, nEigenFuncs);
fprintf('*********************************************************\n');
for r = 1:R
    fprintf([compareTo '. t = %d\n'], r);
    for m = 1:nEigenFuncs
        fprintf('m = %d; Ratio: %f\n', m-1, max(tPhiToCompare(r,:,m))/max(tPhiNumeric(r,:,m)))
        if strcmp(compareTo, 'Analytic')
            mSqNorm(r,m) = sqrt((1/nPoints)*sum((tPhiToCompare(r,:,m).' - tPhiNumeric(r,:,m).').^2));
        elseif strcmp(compareTo, 'Nystrom')
            mSqNorm(r,m) = sqrt((1/nPoints)*sum((tPhiToCompare(r,:,m).' - tPhiNumeric(r,:,m).').^2));
        else
            error('unknown type')
        end
        mSqNorm(r,m) = norm(tPhiToCompare(r,:,m) - tPhiNumeric(r,:,m));
    end
    fprintf('\n');
end
fprintf('*********************************************************\n');
vRMSE = sqrt(sum(mSqNorm,1)/R);
