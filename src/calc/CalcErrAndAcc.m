function [vAcc, vAccStd, vRmse, vMse, vCoh, mErrors] = CalcErrAndAcc(tSig, tSigRef, compareTo)
[n, nSignals, R] = size(tSig);
if isequal(tSig(:,1,1),floor(tSig(:,1,1)))
    errFunc = '0-1';
else
    errFunc = 'norm';
end

mErr = zeros(R, nSignals);
mErrNormed = zeros(R, nSignals);
mCoherence = zeros(R, nSignals);
%fprintf('*********************************************************\n');
for r = 1:R
    %fprintf(['Accuracy of ', compareTo '. r = %d\n'], r);
    mSig = squeeze(tSig(:,:,r));
    mSigRef = squeeze(tSigRef(:,:,r));
    if strcmp(errFunc, 'norm')
        mErr(r,:) = vecnorm(mSig-mSigRef,2);
        mErrNormed(r,:) = min(vecnorm(mSig-mSigRef,2)./vecnorm(mSig,2),1);
        mCoherence(r,:) = (vecnorm(mSig-mSigRef,2).^2)./(vecnorm(mSig,2).*vecnorm(mSigRef,2));
    elseif strcmp(errFunc, '0-1')
        mErr(r,:)       = sum(mSig ~= mSigRef) / n;
        mErrNormed(r,:) = sum(mSig ~= mSigRef) / n;
        mCoherence(r,:) = sum(mSig ~= mSigRef) / n;
    end
    %for m = 1:nSignals
        %fprintf('\t(m=%2d) %6.4f%%\t', m-1, 100*(1-mErrNormed(r,m)))
        %if (mod(m,8) == 0) && m < nSignals
            %fprintf('\n');
        %end
    %end
    %fprintf('\n');
end
% fprintf('*********************************************************\n');

mErrors = mErrNormed.';

% The sum is over the iterations, r=1:R
vRmse = sqrt(sum(mErr.^2,1)/R).';
vRmseStd = std(mErr,[],1).';
vMse = (sum(mErrNormed.^2,1)/R).';
vCoh = (sum(mCoherence,1)/R).';
vAcc = (100*sum((1-mErrNormed),1)/R).';
vAccStd = 100*std(mErrNormed,[],1).';
assert(all(vAcc >= 0) && all(vAccStd >= 0));

%vAcc = vRmse;
%vAccStd = vRmseStd;
end


