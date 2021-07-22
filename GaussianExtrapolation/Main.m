%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')

%% Plot params
% sPlotParams.outputFolder = 'figs';
sPlotParams.b_plotAllEvecs             = false;
sPlotParams.b_GSPBoxPlots              = false;
sPlotParams.b_plotWeights              = true;
sPlotParams.b_plotVRec                 = false;
sPlotParams.b_plotTransDemos           = false;
sPlotParams.b_plotOrigVsInterpEvecs    = true;
sPlotParams.b_plotTildeFiguresForDebug = false;
sPlotParams.b_plotData                 = false;
sPlotParams.b_plotDataVsGmm            = true;
sPlotParams.b_plotInnerProductMatrices = true;
sPlotParams.b_plotHistogram            = false;
sPlotParams.b_plotC                    = true;
%% Number of eigenvectors/eigenfunctions
M                  = 50;
MTilde             = 150; % M
plotInd            = [0,11]; %[M-4, M-1];
%% PolyFit params
b_saturateT        = true;
pCdfDegree         = 10;
invpCdfDegree      = 10;
%% KDE params
kde_bw             = []; %1e-3;
%% GMM params
gmmRegVal          = 1e-3;
gmmMaxIter         = 2000;
gmmNumComponents   = 10; % 1
%% Method parameters
b_debugUseAnalytic = false;
b_applyT           = false;
b_kde              = true;
b_gmmInsteadOfT    = true;
b_forceCtoIdentity = false;
%% Dataset parameters
dim                = 2;
nComponents        = 1;
n                  = 2048;
N                  = 4096;
k                  = 100;
nnValue            = 'ZeroOne'; % 'ZeroOne' / 'Distance'
verticesPDF        = 'MnistLatentVAE'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE'
origGraphAdjacency = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
interpMethod       = 'AddPoints'; % 'NewPoints' / 'AddPoints'
matrixType         = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'NormLap'
%% Verify
assert(~b_debugUseAnalytic || (b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
assert((b_gmmInsteadOfT == ~b_applyT) || (~b_gmmInsteadOfT && ~b_applyT))
%% RMSE
R = 1;
sDatasetParams.xMin = [0 0];
sDatasetParams.xMax = [4 1];

for c = 1:nComponents
    sDatasetParams.mu{c} = 10*(c-1)*ones(1,dim);
    sDatasetParams.sigma{c} = 1*eye(dim);
end

sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, interpMethod, sDatasetParams);
dim            = sDataset.dim;
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
    xTrain   = sDataset.sData.x;
    yTrain   = sDataset.sData.y;
    xInt     = sDataset.sData.xt;
    if r == 1 && sPlotParams.b_plotData && dim <= 3
        PlotDataset(sPlotParams, xTrain, yTrain, verticesPDF, 'Training set');
    end
    % ----------------------------------------------------------------------------------------------
    % Original graph
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [V, lambda] = SimpleCalcAnalyticEigenfunctions(xTrain, omega, sDatasetParams.sigma, sDatasetParams.mu, M);
        lambda = n*lambda;
    else
        [W, dist] = SimpleCalcAdjacency(xTrain, origGraphAdjacency, omega, k, nnValue);
        [V, lambda] = EigsByType(W, M, matrixType);
    end
    if r == 1 && sPlotParams.b_plotWeights
        PlotWeightsMatrix(sPlotParams, W, dist, xTrain, origGraphAdjacency, verticesPDF, omega, k);
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
    if r == 1 && dim <= 2 && sPlotParams.b_plotHistogram
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
        if r == 1 && sPlotParams.b_plotDataVsGmm && dim <= 3
            nGmmPoints = 5000;
            pltTitle = ['Dataset with n = ', num2str(n), ' points'];
            plt2Title = ['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(gmmNumComponents)];
            
            windowStyle = 'normal';
            PlotDataset(sPlotParams, xTrain, yTrain, verticesPDF, pltTitle, sDistParams.GMModel, nGmmPoints, plt2Title, windowStyle);
            
        end
    else
        [PhiTilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
    end
    if b_forceCtoIdentity
        C = zeros(MTilde, M);
        C(1:M,1:M) = eye(M);
    else
        C = pinv(PhiTilde)*V;
    end
    if r == 1 && sPlotParams.b_plotC
        figure('Name', 'C');
        imagesc(C);
        colormap('jet');
        colorbar();
        title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); 
        set(gca,'FontSize', 14);
    end
    
    if r == 1 && sPlotParams.b_plotTildeFiguresForDebug && dim <= 3
        if b_gmmInsteadOfT && ~strcmp(verticesPDF,'Gaussian')
            msgBoxTitle = 'plotTildeFiguresForDebug warning';
            msgBoxMsg = 'gmmInsteadOfT is true and data is not Gaussian, this plot has no meaning...';
            b_waitForOk = true;
            MyMsgBox(msgBoxMsg, msgBoxTitle, b_waitForOk)
        else
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
    % Plot original V for comparison
    % ----------------------------------------------------------------------------------------------
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
%         figTitle = 'Numeric vs. analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
%         PlotSpectrum(sPlotParams, sDataset, [], n*lambdaAnalyticTilde, lambdaNumericTilde, [], ...
%             'n \tilde{\lambda}^{\phi}_m', '\tilde{\lambda}^{v}_m', [], figTitle);
        
        if strcmp(origGraphAdjacency, 'GaussianKernel')
            figTitle = [ 'Eigenvectors of ${\bf W}$ (Gaussian kernel) with $\omega = ' num2str(omega) '\quad n = ' num2str(n) '$'];
        elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
            figTitle = [ 'Eigenvectors of ${\bf W}$ (k-NN) with $k = ' num2str(k) '\quad n = ' num2str(n) '$'];
        end
        figName = 'V';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            V, lambda, '\lambda^{v}', [], figTitle, figName, 'v')
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
    VInt = PhiTildeInt*C;
    VIntToCompare = VInt;
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    normFactor = lambda*sqrt(interpRatio);
    lambdaNys = interpRatio*lambda;
    VNys = B.'*V*diag(1./normFactor);
%     VNys = VNys./norm(VNys(:,1));
%     VNysToCompare = sqrt(interpRatio)*VNys;
    VNysToCompare = VNys;

    
    % ----------------------------------------------------------------------------------------------
    % Calculate reference
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [PhiRef, lambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, sDatasetParams.sigma, sDatasetParams.mu, M);
        VRef = PhiRef;
    else
        WRef = SimpleCalcAdjacency(xInt, origGraphAdjacency, omega, k, nnValue);
        [VRef, lambdaRef] = EigsByType(WRef, M, matrixType);
    end
    VRefToCompare = FlipSign(VInt, VRef);
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VRefToCompare, lambdaRef, '\lambda^{{\bf ref}}', [], ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf ref}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VNysToCompare, lambdaNys, '\lambda^{{\bf nys}}', [], ['Nystrom (N = ', num2str(N), ')'], 'VNys', 'v^{{\bf nys}}');
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
            VIntToCompare, lambdaAnalyticTilde, '\lambda^{{\bf int}}', [], ['Ours (N = ', num2str(N), ')'], 'VInt', 'v^{{\bf int}}');    
    end
    
    % Plot inner product of interpolated eigenvectors
    if r == 1 && sPlotParams.b_plotInnerProductMatrices
        pltTitle = 'VRef - ${\bf V}_{{\bf ref}}^T {\bf V}_{{\bf ref}}$';
        figName = 'VRef';
        PlotInnerProductMatrix(sPlotParams, dim, [], 'IP_Matrix', [], VRef, pltTitle, figName);
        pltTitle = 'VInt - ${\bf V}_{{\bf int}}^T {\bf V}_{{\bf int}}$';
        figName = 'Vint';
        PlotInnerProductMatrix(sPlotParams, dim, [], 'IP_Matrix', [], VInt, pltTitle, figName);
        pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
        figName = 'VNys';
        PlotInnerProductMatrix(sPlotParams, dim, [], 'IP_Matrix', [], VNys, pltTitle, figName);
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
plot((0:M-1)', vRmseInt.', '-o', 'LineWidth', 2, 'DisplayName', 'RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$'), hold on;
plot((0:M-1)', vRmseNys.', '-x', 'LineWidth', 2, 'DisplayName', 'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$')
rmseylim = max([vRmseInt vRmseNys]);
ylim([0 rmseylim]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
set(gca,'FontSize', 14);
% set(0,'DefaultFigureWindowStyle',windowStyle)
