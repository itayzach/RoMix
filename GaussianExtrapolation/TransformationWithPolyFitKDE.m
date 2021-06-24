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
sPlotParams.b_plotTildeFiguresForDebug = true;
sPlotParams.b_plotData                 = true;
sPlotParams.b_plotGmm                  = true;
plotInd                                = [0,3]; %[M-4, M-1];
%% Number of eigenvectors/eigenfunctions
M                  = 50;
MTilde             = 80;
%% PolyFit params
b_saturateT        = true;
pCdfDegree         = 10;
invpCdfDegree      = 10;
%% KDE params
kde_bw             = []; %1e-3;
%% GMM params
gmmRegVal          = 0;
gmmMaxIter         = 1000;
gmmNumComponents   = 16;
%% Method parameters
b_debugUseAnalytic = false;
b_applyT           = false;
b_kde              = true;
b_gmmInsteadOfT    = true;
%% Dataset parameters
dim                = 1;
nComponents        = 1;
n                  = 2000;
N                  = 3000;
k                  = 100;
nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
verticesPDF        = 'TwoMoons'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
origGraphAdjacency = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
interpMethod       = 'AddPoints'; % 'NewPoints' / 'AddPoints'
%% Verify
assert(~b_debugUseAnalytic || (b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
assert((b_gmmInsteadOfT == ~b_applyT) || (~b_gmmInsteadOfT && ~b_applyT))
%% RMSE
R = 1;
for c = 1:nComponents
    sDatasetParams.mu{c} = 10*(c-1)*ones(1,dim);
    sDatasetParams.sigma{c} = 1*eye(dim);
end

sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod, sDatasetParams);
dim            = sDataset.dim;
xTrain         = sDataset.sData.x;
xInt           = sDataset.sData.xt;
xMax           = sDataset.xMax;
xMin           = sDataset.xMin;
n              = length(sDataset.sData.x);
N              = length(sDataset.sData.xt);
omega          = sDataset.recommendedOmega;
omegaTilde     = sDataset.recommendedOmegaTilde;
omegaNys       = sDataset.recommendedOmega;
muTilde        = zeros(1,dim); % was mean(xTrain);
sigmaTilde     = diag(ones(dim,1)); % was diag(std(xTrain));
interpRatio    = N/n;
mVIntToCompare = zeros(R, N, M);
mVNysToCompare = zeros(R, N, M);
mVRefToCompare = zeros(R, N, M);
for r = 1:R
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod, sDatasetParams);
    
    if r == 1 && sPlotParams.b_plotData
        PlotDataset(sPlotParams, xTrain, verticesPDF, 'Training set');
    end
    % ----------------------------------------------------------------------------------------------
    % Original graph
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [V, lambda] = SimpleCalcAnalyticEigenfunctions(xTrain, omega, sDatasetParams.sigma, sDatasetParams.mu, M);
        lambda = n*lambda;
    else
        W = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega, k, nnValue);
        [V, Lambda] = eigs(W,M);
        lambda = diag(Lambda);
    end
    if r == 1 && sPlotParams.b_plotWeights
        PlotWeightsMatrix(sPlotParams, W, xTrain, origGraphAdjacency, verticesPDF, omega, k);
    end
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        if strcmp(origGraphAdjacency, 'GaussianKernel')
            figTitle = [ 'Eigenvectors of ${\bf W}$ (Gaussian kernel) with $\omega = ' num2str(omega) '\quad n = ' num2str(n) '$'];
        elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
            figTitle = [ 'Eigenvectors of ${\bf W}$ (k-NN) with $k = ' num2str(k) '\quad n = ' num2str(n) '$'];
        end
        figName = 'V';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            V, lambda, [], [], figTitle, figName, 'v')
     end
    if ~b_kde && dim <= 3
        % ----------------------------------------------------------------------------------------------
        % Estimate CDF, and model it with polyfit (PolyfitEstCdf)
        % ----------------------------------------------------------------------------------------------
        nEvalPoints = min(700, round(n/10));
        b_plotCdf = true;
        [xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = ...
            PolyfitEstCdf(xTrain, nEvalPoints, pCdfDegree, invpCdfDegree, b_plotCdf);
    else
        pCdf = NaN(dim,pCdfDegree+1);
        invpCdf = NaN(dim,pCdfDegree+1);
        xTrainGrid = xTrain;
        estMarginalCdf_xTrain = NaN(n,dim);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Transform using pCdf (T)
    % ----------------------------------------------------------------------------------------------
    if ~b_applyT
        xTildeTrain = xTrain;
        polyvalCdfMatrix = [];
    else
        [xTildeTrain, polyvalCdfMatrix] = T(pCdf, b_saturateT, muTilde, sigmaTilde, xTrain, xTrain, b_kde, kde_bw);
    end
    if r == 1 && dim <= 2
        PlotHistogram(sPlotParams, xTildeTrain, 'Gaussian', 'Histogram of $\tilde{X}$', false);
    end
    if r == 1 && strcmp(verticesPDF, 'Gaussian')
        PlotGaussianSanity(xTrain, xTildeTrain, muTilde, sigmaTilde, ...
            ~b_applyT, xTrainGrid, estMarginalCdf_xTrain, polyvalCdfMatrix, b_kde);
    end

    if r == 1 && sPlotParams.b_plotTransDemos && b_applyT
        for d = 1:dim
            PlotPolyCdfDemonstration1(xMin(d), xMax(d), pCdf(d,:), xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), muTilde(d), sigmaTilde(d,d), b_kde, kde_bw);
            PlotPolyCdfDemonstration2(xMin(d), xMax(d), pCdf(d,:), invpCdf(d,:), muTilde(d), sigmaTilde(d,d), xTrain, b_kde, kde_bw);
        end
        if dim == 2
            PlotPolyCdfDemonstration3_2D(xTrain,xTildeTrain)
        end
    end
    % ----------------------------------------------------------------------------------------------
    % Learn eigenvectors to eigenfunctions transformation (C)
    % ----------------------------------------------------------------------------------------------
    if b_gmmInsteadOfT
        sDistParams = EstimateDistributionParameters(xTildeTrain, gmmNumComponents, gmmRegVal, gmmMaxIter);
        sKernelParams = GetKernelParams(sDistParams, omegaTilde);
        [sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
            = CalcAnalyticEigenvalues(MTilde, sKernelParams, dim, gmmNumComponents);
        [ PhiTilde, lambdaAnalyticTilde ] = ...
            CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeTrain, true);
        if sPlotParams.b_plotGmm
            nGmmPoints = 5000;
            pltTitle = ['Dataset with n = ', num2str(n), ' points'];
            plt2Title = ['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(gmmNumComponents)];
            PlotDataset(sPlotParams, xTrain, verticesPDF, pltTitle, sDistParams.GMModel, nGmmPoints, plt2Title);
        end
    else
        [PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
    end
    C = pinv(PhiTilde)*V;
    if r == 1
        figure('Name', 'C');
        imagesc(C);
        colorbar();
        title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
        set(gca,'FontSize', 14);
    end
    
    if r == 1 && sPlotParams.b_plotTildeFiguresForDebug && dim <= 3
        sqrtnLambdaPhiTilde = sqrt(n*lambdaAnalyticTilde)'.*PhiTilde;
        % Build W tilde and numeric eigenvectors for comparison
        distTilde = pdist2(xTildeTrain, xTildeTrain);
        WTilde = exp(-distTilde.^2/(2*omegaTilde^2));
        [VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
        lambdaNumericTilde = diag(LambdaNumericTilde);
        sqrtLambdaVTilde = abs(sqrt(lambdaNumericTilde))'.*VTilde;
        VTilde = FlipSign(PhiTilde, VTilde);

        figTitle = 'Eigenfunctions vs. eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        figName = 'PhiTilde_VTilde';
        PlotEigenfuncvecScatter(sPlotParams, 'Gaussian', xTildeTrain, [], plotInd(1), plotInd(end), ...
            sqrtnLambdaPhiTilde, [], [], [], figTitle, figName, ...
            '\sqrt{n \tilde{\lambda}^{\phi}}\tilde{\phi}', sqrtLambdaVTilde, '\sqrt{\tilde{\lambda}^{v}}\tilde{v}')
        figTitle = 'Numeric vs. analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        PlotSpectrum(sPlotParams, sDataset, [], n*lambdaAnalyticTilde, lambdaNumericTilde, [], ...
            'n \tilde{\lambda}^{\phi}_m', '\tilde{\lambda}^{v}_m', [], figTitle);
        vPrTilde = SimpleEstPorbablityArea(xTildeTrain, sigmaTilde, muTilde);
        pltTitle = 'Analytic - $\int \phi_i(x) \phi_j(x) p(x) dx = n^d \Phi^T$diag(Pr)$\Phi$';
        figName = 'PhiTilde';
        PlotInnerProductMatrix(sPlotParams, dim, vPrTilde, 'IP_Matrix', [], PhiTilde, pltTitle, figName);
    end
    if r == 1 && sPlotParams.b_plotVRec && dim <= 3
        VRec = PhiTilde*C;

        figTitle = '${\bf V}$ in terms of ${\bf \tilde{\Phi}} \quad ({\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C})$';
        figName = 'VRec';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            VRec, lambda, [], [], figTitle, figName, 'v^{{\bf rec}}')

        b_plotErrVsNodeInd = true;
        PlotEigenDiffs(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), VRec, V, ...
            'VRec_vs_V', 'v^{{\bf rec}}', 'v', b_plotErrVsNodeInd)
        
    end
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with our method
    % ----------------------------------------------------------------------------------------------
    if ~b_applyT
        xTildeInt = xInt;
    else
        xTildeInt = T(pCdf, b_saturateT, muTilde, sigmaTilde, xInt, xTrain, b_kde, kde_bw);
    end
    
    if b_gmmInsteadOfT
        [PhiTildeInt, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeInt, true);
    else
        [PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);
    end
    VInt = sqrt(interpRatio)*PhiTildeInt*C;
    VIntToCompare = abs(VInt);
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    VNys = B.'*V*diag(1./lambda);
    VNysToCompare = abs(VNys);
    
    % ----------------------------------------------------------------------------------------------
    % Calculate reference
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [PhiRef, lambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, sDatasetParams.sigma, sDatasetParams.mu, M);
        VRef = PhiRef;
    else
        WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega, k, nnValue);
        [VRef, LambdaRef] = eigs(WRef,M);
        lambdaRef = diag(LambdaRef);
    end
    VRef = FlipSign(VInt, VRef);
    VRefToCompare = sqrt(interpRatio)*abs(VRef);
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VIntToCompare, [], [], [], ['Ours vs. Reference (N = ', num2str(N), ')'], 'VInt_vs_VRef', 'v^{{\bf int}}', VRefToCompare, 'v^{{\bf ref}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VNysToCompare, [], [], [], ['Nystrom vs. Reference (N = ', num2str(N), ')'], 'VNys_vs_VRef', 'v^{{\bf nys}}', VRefToCompare, 'v^{{\bf ref}}');
    end
    
    % Plot inner product of interpolated eigenvectors
    pltTitle = 'VInt - ${\bf V}_{{\bf int}}^T {\bf V}_{{\bf int}}$';
    figName = 'Vint';
    PlotInnerProductMatrix(sPlotParams, dim, [], 'IP_Matrix', [], VInt/sqrt(interpRatio), pltTitle, figName);
    pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
    figName = 'VNys';
    PlotInnerProductMatrix(sPlotParams, dim, [], 'IP_Matrix', [], VNys/sqrt(interpRatio), pltTitle, figName);
    
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
plot((0:M-1)', vRmseInt.', '-o', 'LineWidth', 2, 'DisplayName', 'RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$'), hold on;
plot((0:M-1)', vRmseNys.', '-x', 'LineWidth', 2, 'DisplayName', 'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$')
rmseylim = max([vRmseInt vRmseNys]);
ylim([0 rmseylim]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);
% set(0,'DefaultFigureWindowStyle',windowStyle)
