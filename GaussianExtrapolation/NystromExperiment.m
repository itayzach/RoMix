%% Restart run
clear; 
close all; 
clc;
rng('default'); % For reproducibility
%% Parameters
vNysRatio = [0.05, 0.4, 0.8];
sParams = GetParameters();

%% Run
T = 50;
for simIdx = 1:2
    if simIdx == 1
        b_useTwoMoonsDataset = false;
    else
        b_useTwoMoonsDataset = true;
    end
    for nysRatio = vNysRatio
        mAnaVsNum = zeros(T, sParams.PlotSpectM);
        mNysVsNum = zeros(T, sParams.PlotSpectM);
        for t = 1:T
            sParams = GetParameters(b_useTwoMoonsDataset, nysRatio);
            [ mPhiAnalytic, vLambdaAnalytic ]  = CalcAnalyticEigenfunctions(sParams, sParams.x_rand);
            [ mPhiNumeric, vLambdaNumeric ]    = CalcNumericEigenvectors(sParams);
            [ mPhiNystrom, vLambdaNystrom ]    = CalcNystromEigenvectors(sParams);
            
            mAnaVsNum(t,:) = vecnorm(mPhiAnalytic - mPhiNumeric).^2;
            mNysVsNum(t,:) = vecnorm(mPhiNystrom - mPhiNumeric).^2;
        end
        vAnaVsNum = sqrt(sum(mAnaVsNum,1)/T);
        vNysVsNum = sqrt(sum(mNysVsNum,1)/T);
        PlotEigenDiffs(sParams, vAnaVsNum, vNysVsNum)
        
        % Plot histogram and eigenvectors of the last iteration
        PlotHistogram(sParams);
        PlotEigenfuncvecScatter(sParams, 0, 19, mPhiAnalytic, vLambdaAnalytic, 'Analytic')
        PlotEigenfuncvecScatter(sParams, 0, 19, mPhiNumeric, vLambdaNumeric, 'Numeric')
        PlotEigenfuncvecScatter(sParams, 0, 19, mPhiNystrom, vLambdaNystrom, 'Nystrom')
    end
end