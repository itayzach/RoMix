function sClassifier = BuildClassifier(sSimParams, sDataset, sKernelParams, name, gamma_A, gamma_I)
sClassifier.name = name;
sClassifier.gamma_A = gamma_A;
sClassifier.gamma_I = gamma_I;
sClassifier.omega = sKernelParams.omega;

b_normalize = false;

%% Train
[ mPhiAnalyticTrain, vLambdaAnalytic ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.x, b_normalize);

W = CalcAdjacency(sKernelParams, sDataset.sData.x);
D = diag(sum(W,1));
L = D - W;

c = eigrls(sDataset.sData.y, mPhiAnalyticTrain, diag(vLambdaAnalytic), gamma_A, gamma_I, L);
sClassifier.c = c;
sClassifier.vPhi_xTrain_c = mPhiAnalyticTrain*c;

%% Test
[ mPhiAnalyticTest, ~ ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.sData.xt, b_normalize);
sClassifier.vPhi_xTest_c = mPhiAnalyticTest*c;
vPredictions = sign(sClassifier.vPhi_xTest_c);
vTestIdx = find(sDataset.sData.yt);
sClassifier.error = 100*sum(vPredictions(vTestIdx) ~= sDataset.sData.yt(vTestIdx)) / length(vTestIdx);
fprintf('error: %f\n', sClassifier.error);
%% Extrapolate
xMax = max(max(sDataset.sData.x,[],1));
xMin = min(min(sDataset.sData.x,[],1));
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[sClassifier.XX1,sClassifier.XX2] = meshgrid(x1,x2);
X=[sClassifier.XX1(:) sClassifier.XX2(:)];    

tic;
[ mPhi_m_X, ~ ] = CalcAnalyticEigenfunctions(sSimParams.CalcEigenFuncsM, sKernelParams, X, b_normalize);
run_time = toc;
fprintf('Phi time = %f\n', run_time);

if sSimParams.b_plotEigenfunctions
    firstEigenIdxToPlot = 0;
    lastEigIdxToPlot = sSimParams.PlotEigenFuncsM-1;
    PlotEigenfuncvecScatter(sSimParams, 'Two_moons', X, [], firstEigenIdxToPlot, lastEigIdxToPlot, mPhi_m_X, vLambdaAnalytic, 'Analytic',sClassifier.c)
end

sClassifier.mPhi_X_c = reshape(mPhi_m_X*c, length(x1), length(x2));

end