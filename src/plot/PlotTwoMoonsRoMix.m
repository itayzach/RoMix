function PlotTwoMoonsRoMix(sPlotParams, sDataset, sKernelParams, mSigRecPhi, mSigHatPhi, gamma1, gamma2)
MTilde = size(mSigHatPhi,1);
sClassifier.MTilde = MTilde;
sClassifier.omega = sKernelParams.omega;
sClassifier.gamma_A = gamma1;
sClassifier.gamma_I = gamma2;
sClassifier.vPhi_xTrain_c = mSigRecPhi;
xMax = max([sDataset.sData.x(:); sDataset.sData.xt(:)]);
xMin = min([sDataset.sData.x(:); sDataset.sData.xt(:)]);
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[sClassifier.XX1,sClassifier.XX2] = meshgrid(x1,x2);
xMeshGrid = [sClassifier.XX1(:) sClassifier.XX2(:)];
PhiMeshgrid = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xMeshGrid);
sClassifier.mPhi_X_c = reshape(PhiMeshgrid*mSigHatPhi, length(x1), length(x2));
PlotClassifier(sPlotParams, sDataset, sClassifier);
end