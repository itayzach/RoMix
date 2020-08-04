function vRMSE = CalcRMSE(sSimParams, T, tPhi1, tPhi2)
mSqNorm = zeros(T, sSimParams.PlotEigenFuncsM);

for t = 1:T
    mSqNorm(t,:) = vecnorm(tPhi1(t,:,:) - tPhi2(t,:,:)).^2;
end
vRMSE = sqrt(sum(mSqNorm,1)/T);