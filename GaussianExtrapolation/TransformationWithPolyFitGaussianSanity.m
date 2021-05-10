%% Restart
clc; clear; rng('default');
close all;
set(0,'DefaultFigureWindowStyle','normal')
sSimParams.outputFolder = 'figs';
%% Random data
dim = 1;
nComponents = 1;
n = 1000;
N = 5000;
M = 30;
MTilde = 30;
b_saturateT = true;
pCdfDegree = 10;
invpCdfDegree = 10;
verticesPDF = 'Gaussian'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
origGraphAdjacency = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
interpMethod = 'AddPoints'; % 'NewPoints' / 'AddPoints'

%% Debug
b_debugUseAnalytic = false;
b_debugUseNormalCDF = true;
%% Params from sDataset
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod);

%% RMSE
R = 1;
mVIntToCompare = zeros(R, N, M);
mVNysToCompare = zeros(R, N, M);
mVRefToCompare = zeros(R, N, M);
for r = 1:R
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod);
    xTrain = sDataset.sData.x;
    xInt = sDataset.sData.xt;
    n = length(sDataset.sData.x);
    N = length(sDataset.sData.xt);
    omega = sDataset.recommendedOmega;
    omegaTilde = omega;
    muTilde = zeros(dim,1); % was mean(xTrain);
    sigmaTilde = diag(ones(dim,1)); % was diag(std(xTrain));
    interpRatio = N/n;
    
    % ----------------------------------------------------------------------------------------------
    % Original graph
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        mu = 0; sigma = 1;
        [V, lambda] = SimpleCalcAnalyticEigenfunctions(xTrain, omega, sigma, mu, M);
        lambda = n*lambda;
    else
        W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega);
        [V, Lambda] = eigs(W,M);
        lambda = diag(Lambda);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Learn pCdf, calc eigen functions and C
    % ----------------------------------------------------------------------------------------------
    bins = min(700, n);
    [xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = PolyfitEstCdf(xTrain, bins, pCdfDegree, invpCdfDegree);
    [xTildeTrain, polyvalCdf] = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain, b_debugUseNormalCDF, xTrainGrid, estMarginalCdf_xTrain);
    
%     PlotPolyCdfDemonstration1(sDataset.xMin, sDataset.xMax, pCdf, xTrainGrid, estMarginalCdf_xTrain, muTilde, sigmaTilde);
%     PlotPolyCdfDemonstration2(sDataset.xMin, sDataset.xMax, pCdf, invpCdf, muTilde, sigmaTilde);
    
    [PhiTilde, ~] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
    C = pinv(PhiTilde)*V;
    if r == 1
        figure('Name', 'C');
        imagesc(C);
        colorbar();
        title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
        set(gca,'FontSize', 14);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with our method
    % ----------------------------------------------------------------------------------------------
    xTildeInt = T(pCdf, true, muTilde, sigmaTilde, xInt, b_debugUseNormalCDF, xTrainGrid, estMarginalCdf_xTrain);
    [PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);
    VInt = sqrt(interpRatio)*PhiTildeInt*C;
    VIntToCompare = abs(VInt);
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omega^2));
    VNys = B.'*V*diag(1./lambda);
    VNysToCompare = abs(VNys);
    
    % ----------------------------------------------------------------------------------------------
    % Calculate reference
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [PhiRef, ~] = SimpleCalcAnalyticEigenfunctions(xInt, omega, sigma, mu, M);
        VRef = PhiRef;
    else
        WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega);
        [VRef, ~] = eigs(WRef,M);
    end
    VRef = FlipSign(VInt, VRef);
    VRefToCompare = sqrt(interpRatio)*abs(VRef);
    if r == 1
        figTitle = 'Ours vs. Reference';
        figName = 'VInt_vs_VRef';
        PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], 0, 3, ...
            VIntToCompare, [], [], [], figTitle, figName, 'v^{{\bf int}}', VRefToCompare, 'v^{{\bf ref}}');
        
        figTitle = 'Nystrom vs. Reference';
        figName = 'VNys_vs_VRef';
        PlotEigenfuncvecScatter(sSimParams, verticesPDF, xInt, [], 0, 3, ...
            VNysToCompare, [], [], [], figTitle, figName, 'v^{{\bf nys}}', VRefToCompare, 'v^{{\bf ref}}');
    end
    
    % Save
    mVIntToCompare(r,:,:) = VIntToCompare;
    mVNysToCompare(r,:,:) = VNysToCompare;
    mVRefToCompare(r,:,:) = VRefToCompare;
    
end
vRmseInt = CalcRMSE(mVIntToCompare, mVRefToCompare, 'Analytic');
vRmseNys = CalcRMSE(mVNysToCompare, mVRefToCompare, 'Nystrom');

windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')
figure;
plot((0:M-1)', [vRmseInt.' vRmseNys.'], 'LineWidth', 2)
ylim([0 0.5]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Ours', 'Nystrom', 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);
set(0,'DefaultFigureWindowStyle',windowStyle)
