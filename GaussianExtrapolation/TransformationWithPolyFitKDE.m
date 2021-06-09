%% Restart
clc; clear; rng('default');
close all;
set(0,'DefaultFigureWindowStyle','docked')

%% Plot params
sPlotParams.outputFolder = 'figs';
sPlotParams.b_plotAllEvecs             = false;
sPlotParams.b_GSPBoxPlots              = false;
sPlotParams.b_plotLaplacianEvecs       = false;
sPlotParams.b_plotWeights              = false;
sPlotParams.b_plotVRec                 = false;
sPlotParams.b_plotTransDemos           = true;
sPlotParams.b_plotOrigVsInterpEvecs    = true;
plotInd                                = [0,3]; %[M-4, M-1];
%% Parameters
% Method parameters
M                  = 20;
MTilde             = 50;
b_saturateT        = true;
pCdfDegree         = 10;
invpCdfDegree      = 10;
b_kde              = true;
% Dataset parameters
dim                = 1;
nComponents        = 1;
n                  = 5000;
N                  = 8000;
verticesPDF        = 'Grid'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
origGraphAdjacency = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
interpMethod       = 'AddPoints'; % 'NewPoints' / 'AddPoints'
%% Debug
sDebugParams.b_debugUseAnalytic = false;
sDebugParams.b_debugUseNormalCDF = false;

assert(~sDebugParams.b_debugUseAnalytic || (sDebugParams.b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
assert(~sDebugParams.b_debugUseNormalCDF || (sDebugParams.b_debugUseNormalCDF && strcmp(verticesPDF,'Gaussian')))
%% RMSE
R = 1;
sDatasetParams.mu = 0*ones(1,dim);
sDatasetParams.sigma = 1*eye(dim);
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod, sDatasetParams);
N = length(sDataset.sData.xt);
mVIntToCompare = zeros(R, N, M);
mVNysToCompare = zeros(R, N, M);
mVRefToCompare = zeros(R, N, M);
for r = 1:R
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    sDatasetParams.mu = 0*ones(1,dim);
    sDatasetParams.sigma = 1*eye(dim);
    sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod, sDatasetParams);
    dim            = sDataset.dim;
    xTrain         = sDataset.sData.x;
    xInt           = sDataset.sData.xt;
    xMax           = sDataset.xMax;
    xMin           = sDataset.xMin;
    n              = length(sDataset.sData.x);
    N              = length(sDataset.sData.xt);
    sDatasetParams = sDataset.sDatasetParams;
    omega          = sDataset.recommendedOmega;
    omegaTilde     = sDataset.recommendedOmega;
    muTilde        = zeros(1,dim); % was mean(xTrain);
    sigmaTilde     = diag(ones(dim,1)); % was diag(std(xTrain));
    interpRatio    = N/n;
    % ----------------------------------------------------------------------------------------------
    % Original graph
    % ----------------------------------------------------------------------------------------------
    if sDebugParams.b_debugUseAnalytic
        [V, lambda] = SimpleCalcAnalyticEigenfunctions(xTrain, omega, sDatasetParams.sigma, sDatasetParams.mu, M);
        lambda = n*lambda;
    else
        W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega);
        [V, Lambda] = eigs(W,M);
        lambda = diag(Lambda);
    end
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs)
        if strcmp(origGraphAdjacency, 'GaussianKernel')
        figTitle = [ 'Eigenvectors of ${\bf W}$ (Gaussian kernel) with $\omega = ' num2str(omega) '\quad n = ' num2str(n) '$'];
        elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
            figTitle = [ 'Eigenvectors of ${\bf W}$ (k-NN) with $k = ' num2str(k) '\quad n = ' num2str(n) '$'];
        end
        figName = 'V';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            V, lambda, [], [], figTitle, figName, 'v')
    end
    if ~b_kde
        % ----------------------------------------------------------------------------------------------
        % Estimate CDF, and model it with polyfit (PolyfitEstCdf)
        % ----------------------------------------------------------------------------------------------
        nEvalPoints = min(700, round(n/10));
        b_plotCdf = true;
        [xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = ...
            PolyfitEstCdf(xTrain, nEvalPoints, pCdfDegree, invpCdfDegree, false, b_plotCdf);
    else
        pCdf = [];
        invpCdf = []; 
        xTrainGrid = xTrain;
        for d=1:dim
            estMarginalCdf_xTrain(:,d) = ksdensity(xTrain(:,d), xTrain(:,d), 'Function', 'cdf');
        end
    end
    
    % ----------------------------------------------------------------------------------------------
    % Transform using pCdf (T)
    % ----------------------------------------------------------------------------------------------
    if sDebugParams.b_debugUseNormalCDF
        xTildeTrain = xTrain;
        polyvalCdfMatrix = [];
    else
        [xTildeTrain, polyvalCdfMatrix] = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain, xTrain, b_kde);
    end
    if r == 1 && dim <= 2
        PlotHistogram(sPlotParams, xTildeTrain, 'Gaussian', 'Histogram of $\tilde{X}$', false);
    end
    if r == 1 && strcmp(verticesPDF, 'Gaussian')
        PlotGaussianSanity(xTrain, xTildeTrain, muTilde, sigmaTilde, ...
            sDebugParams.b_debugUseNormalCDF, xTrainGrid, estMarginalCdf_xTrain, polyvalCdfMatrix);
    end

    if r == 1 && sPlotParams.b_plotTransDemos && ~b_kde
        for d = 1:dim
            PlotPolyCdfDemonstration1(xMin(d), xMax(d), pCdf(d,:), xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), muTilde(d), sigmaTilde(d,d));
%             PlotPolyCdfDemonstration2(xMin(d), xMax(d), pCdf(d,:), invpCdf(d,:), muTilde(d), sigmaTilde(d,d));
        end
        if dim == 2
            PlotPolyCdfDemonstration3_2D(xTrain,xTildeTrain)
        end
    end
    % ----------------------------------------------------------------------------------------------
    % Learn eigenvectors to eigenfunctions transformation (C)
    % ----------------------------------------------------------------------------------------------
    [PhiTilde, ~] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
    C = pinv(PhiTilde)*V;
    if r == 1
        figure('Name', 'C');
        imagesc(C);
        colorbar();
        title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
        set(gca,'FontSize', 14);
    end
    
    if r == 1 && sPlotParams.b_plotVRec
        VRec = PhiTilde*C;

        figTitle = '${\bf V}$ in terms of ${\bf \tilde{\Phi}} \quad ({\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C})$';
        figName = 'VRec';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            VRec, lambda, [], [], figTitle, figName, 'v^{{\bf rec}}')

        b_plotErrVsNodeInd = true;
        PlotEigenDiffs(sPlotParams, sDataset, [], plotInd(1), plotInd(end), VRec, V, 'VRec_vs_V', 'v^{{\bf rec}}', 'v', b_plotErrVsNodeInd)
    end
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with our method
    % ----------------------------------------------------------------------------------------------
    if sDebugParams.b_debugUseNormalCDF
        xTildeInt = xInt;
    else
        xTildeInt = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain, xInt, true);
    end  
    
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
    if sDebugParams.b_debugUseAnalytic
        [PhiRef, lambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, sDatasetParams.sigma, sDatasetParams.mu, M);
        VRef = PhiRef;
    else
        WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega);
        [VRef, LambdaRef] = eigs(WRef,M);
        lambdaRef = diag(LambdaRef);
    end
    VRef = FlipSign(VInt, VRef);
    VRefToCompare = sqrt(interpRatio)*abs(VRef);
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs)
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VIntToCompare, [], [], [], 'Ours vs. Reference', 'VInt_vs_VRef', 'v^{{\bf int}}', VRefToCompare, 'v^{{\bf ref}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VNysToCompare, [], [], [], 'Nystrom vs. Reference', 'VNys_vs_VRef', 'v^{{\bf nys}}', VRefToCompare, 'v^{{\bf ref}}');
%         PlotEigenDiffs(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(end), VIntToCompare, VRefToCompare, ...
%             'VInt_vs_VRef', 'v^{{\bf int}}', 'v^{{\bf ref}}', true)
%         PlotEigenDiffs(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(end), VNysToCompare, VRefToCompare, ...
%             'VNys_vs_VRef', 'v^{{\bf nys}}', 'v^{{\bf ref}}', true)
    end
    
    % Save
    mVIntToCompare(r,:,:) = VIntToCompare;
    mVNysToCompare(r,:,:) = VNysToCompare;
    mVRefToCompare(r,:,:) = VRefToCompare;
    
end
vRmseInt = CalcRMSE(mVIntToCompare, mVRefToCompare, 'Analytic');
vRmseNys = CalcRMSE(mVNysToCompare, mVRefToCompare, 'Nystrom');

% windowStyle = get(0,'DefaultFigureWindowStyle');
% set(0,'DefaultFigureWindowStyle','normal')
figure('Name', 'RMSE');
plot((0:M-1)', [vRmseInt.' vRmseNys.'], 'LineWidth', 2)
rmseylim = max([vRmseInt vRmseNys]);
ylim([0 rmseylim]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', ...
    'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
    'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);
% set(0,'DefaultFigureWindowStyle',windowStyle)
