function vRMSE = CalcRMSE(sSimParams, tPhi1, tPhi2)

T = size(tPhi1, 1);

mSqNorm = zeros(T, sSimParams.CalcEigenFuncsM);

for t = 1:T
    mSqNorm(t,:) = (vecnorm(tPhi1(t,:,:) - tPhi2(t,:,:))./vecnorm(tPhi2(t,:,:))).^2;
end
vRMSE = sqrt(sum(mSqNorm,1)/T);