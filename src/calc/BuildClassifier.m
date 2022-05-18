function sClassifier = BuildClassifier(sSimParams, sDataset, sKernelParams, name, gamma_A, gamma_I)
sClassifier.name = name;
sClassifier.gamma_A = gamma_A;
sClassifier.gamma_I = gamma_I;
sClassifier.omega = sKernelParams.omega;
sClassifier.MTilde = sSimParams.CalcEigenFuncsM;

%% Train
[ mPhiAnalyticTrain, vLambdaAnalytic ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.x);

dist = pdist2(sDataset.sData.x, sDataset.sData.x);
W = exp(-dist.^2/(2*sKernelParams.omega^2));

d = sum(W,1);
I = eye(length(sDataset.sData.x));
Ln = I - diag(d.^(-0.5))*W*diag(d.^(-0.5));

PlotWeightsMatrix(sSimParams, W, [], diag(d), sDataset.sData.x, 'GaussianKernel', sKernelParams.omega);

c = eigrls(sDataset.sData.y, mPhiAnalyticTrain, diag(vLambdaAnalytic), gamma_A, gamma_I, Ln);
sClassifier.c = c;
sClassifier.vPhi_xTrain_c = mPhiAnalyticTrain*c;

%% Test
[ mPhiAnalyticTest, ~ ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.xt);
sClassifier.vPhi_xTest_c = mPhiAnalyticTest*c;
vPredictions = sign(sClassifier.vPhi_xTest_c);
vTestIdx = find(sDataset.sData.yt);
sClassifier.error = 100*sum(vPredictions(vTestIdx) ~= sDataset.sData.yt(vTestIdx)) / length(vTestIdx);
fprintf('error: %f\n', sClassifier.error);
%% Extrapolate
xMax = max([sDataset.sData.x(:); sDataset.sData.xt(:)]);
xMin = min([sDataset.sData.x(:); sDataset.sData.xt(:)]);
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[sClassifier.XX1,sClassifier.XX2] = meshgrid(x1,x2);
X=[sClassifier.XX1(:) sClassifier.XX2(:)];    

tic;
[ mPhi_m_X, ~ ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, X);
run_time = toc;
fprintf('Phi time = %.2f ms\n', run_time*1e3);

if sSimParams.b_plotEigenfunctions
    firstEigenIdxToPlot = 0;
    lastEigIdxToPlot = sSimParams.PlotEigenFuncsM-1;
    figTitle = 'Eigenfunctions';
    figName = 'Eigenfunctions';
    PlotEigenfuncvecScatter(sSimParams, 'Two_moons', X, [], firstEigenIdxToPlot, lastEigIdxToPlot, ...
            mPhi_m_X, vLambdaAnalytic, '\lambda', [], figTitle, figName, '\phi');
end

sClassifier.mPhi_X_c = reshape(mPhi_m_X*c, length(x1), length(x2));
sClassifier.vLambdaAnalytic = vLambdaAnalytic;
end