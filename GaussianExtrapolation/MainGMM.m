%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')
%% Plot params
sPlotParams = GetPlotParams();
%% Dataset parameters
dim                 = 2;
n                   = 2048;
N                   = 4096;
k                   = round(0.01*N);
nnValue             = 'ZeroOne'; % 'ZeroOne' / 'Distance'
verticesPDF         = 'MnistLatentVAE'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE'
adjacencyType       = 'GaussianKernel'; % 'NearestNeighbor' / 'GaussianKernel'
matrixForEigs       = 'Adjacency'; % 'Adjacency' / 'RandomWalk' / 'Laplacian' / 'NormLap'
%% Params for Grid/Uniform
sDatasetParams.xMin = [0 0];
sDatasetParams.xMax = [4 1];
%% Params for Gaussian
nComponents = 1;
for c = 1:nComponents
    sDatasetParams.mu{c} = 10*(c-1)*ones(1,dim);
    sDatasetParams.sigma{c} = 1*eye(dim);
end
%% Number of eigenvectors/eigenfunctions
M                  = 50;
MTilde             = 500; % M
%% GMM params
gmmRegVal          = 1e-3;
gmmMaxIter         = 2000;
gmmNumComponents   = 16; % 1
%% Method parameters
b_debugUseAnalytic = false;
b_forceCtoIdentity = false;
%% Verify
assert(~b_debugUseAnalytic || (b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
%% Dataset parameters
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, 'AddPoints', sDatasetParams);
dim            = sDataset.dim;
xMax           = sDataset.xMax;
xMin           = sDataset.xMin;
n              = length(sDataset.sData.x);
N              = length(sDataset.sData.xt);
interpRatio    = N/n;
omega          = sDataset.recommendedOmega;
omegaTilde     = sDataset.recommendedOmegaTilde;
omegaNys       = sDataset.recommendedOmega;
%% plotInd
if dim == 1
    plotInd = [0,4];
else
    plotInd = [0,11];
end
%% Run
R = 1;
mVIntToCompare = zeros(R, N, M);
mVNysToCompare = zeros(R, N, M);
mVRefToCompare = zeros(R, N, M);
for r = 1:R
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, 'AddPoints', sDatasetParams);
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
        [V, adjLambda] = SimpleCalcAnalyticEigenfunctions(xTrain, omega, ...
            sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        adjLambda = n*adjLambda;
    else
        [W, dist] = SimpleCalcAdjacency(xTrain, adjacencyType, omega, k, nnValue);
        [V, adjLambda, matLambda] = EigsByType(W, M, matrixForEigs);
        if r == 1 && sPlotParams.b_plotWeights
            PlotWeightsMatrix(sPlotParams, W, dist, xTrain, adjacencyType, verticesPDF, omega, k);
        end
    end

    % ----------------------------------------------------------------------------------------------
    % Calculate lambdaAnalyticTilde and PhiTilde(xTildeTrain)
    % ----------------------------------------------------------------------------------------------
    xTildeTrain = xTrain;
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
    % ----------------------------------------------------------------------------------------------
    % Learn eigenvectors to eigenfunctions transformation (C)
    % ----------------------------------------------------------------------------------------------
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
    
    % ----------------------------------------------------------------------------------------------
    % Plot PhiTilde
    % ----------------------------------------------------------------------------------------------
    if r == 1 && sPlotParams.b_plotTildeFiguresForDebug && dim <= 3
        if ~strcmp(verticesPDF,'Gaussian')
            figTitle = 'Eigenfunctions of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
            figName = 'PhiTilde_VTilde';
            PlotEigenfuncvecScatter(sPlotParams, 'Gaussian', xTildeTrain, [], 40, 50, ...
                PhiTilde, lambdaAnalyticTilde, '\tilde{\lambda}^{{\phi}}', [], ...
                figTitle, figName, '\tilde{\phi}' )
            figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
            PlotSpectrum(sPlotParams, sDataset, [], lambdaAnalyticTilde, [], [], ...
                '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
%             pltTitle = 'Analytic - $\int \phi_i(x) \phi_j(x) p(x) dx = n^d \Phi^T$diag(Pr)$\Phi$';
%             figName = 'PhiTilde';
%             PlotInnerProductMatrix(sPlotParams, dim, vPrTilde, 'IP_Matrix', [], PhiTilde(:,1:10), pltTitle, figName);
        end
    end
    
    % ----------------------------------------------------------------------------------------------
    % Plot original V for comparison
    % ----------------------------------------------------------------------------------------------
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        figTitle = ['Eigenvectors of ', matrixForEigs];
        if strcmp(adjacencyType, 'GaussianKernel')
            figTitle2 = [' (generated from Gaussian kernel with $\omega = ', num2str(omega), '$)'];
        elseif strcmp(adjacencyType, 'NearestNeighbor')
            figTitle2 = [' (generated from k-NN with $k = ', num2str(k), '$)'];
        end
        figTitle3 = [' $\quad n = ', num2str(n), '$ points'];
        figName = 'V';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            V, matLambda, '\lambda^{v}', [], [figTitle, figTitle2, figTitle3], figName, 'v')
     end

    % ----------------------------------------------------------------------------------------------
    % Interpolate with our method
    % ----------------------------------------------------------------------------------------------
    xTildeInt = xInt;
    [PhiTildeInt, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeInt, true);
    VInt = PhiTildeInt*C;
    VIntToCompare = VInt;
    
    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = pdist2(xTrain, xInt);
    B = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    normFactor = adjLambda*sqrt(interpRatio);
    VNys = B.'*V*diag(1./normFactor);
    if ~strcmp(matrixForEigs, 'Adjacency')
        VNys = VNys./norm(VNys(:,1));
    end
    VNysToCompare = VNys;
    lambdaNys = interpRatio*adjLambda;

    
    % ----------------------------------------------------------------------------------------------
    % Calculate reference
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [PhiRef, adjLambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, ...
            sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        VRef = PhiRef;
    else
        WRef = SimpleCalcAdjacency(xInt, adjacencyType, omega, k, nnValue);
        [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, M, matrixForEigs);
    end
    VRefToCompare = FlipSign(VInt, VRef);
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        if dim == 1
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VRefToCompare, [], [], [], ['Reference vs. Nystrom (N = ', num2str(N), ')'], ...
                'VRef', 'v^{{\bf ref}}', VNysToCompare, 'v^{{\bf nys}}');
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VRefToCompare, [], [], [], ['Reference vs. Ours (N = ', num2str(N), ')'], ...
                'VRef', 'v^{{\bf ref}}', VIntToCompare, 'v^{{\bf int}}');
        else
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VRefToCompare, matLambdaRef, '\lambda^{{\bf ref}}', [], ...
                ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf ref}}');
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VNysToCompare, lambdaNys, '\lambda^{{\bf nys}}', [], ...
                ['Nystrom (N = ', num2str(N), ')'], 'VNys', 'v^{{\bf nys}}');
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VIntToCompare, lambdaAnalyticTilde, '\lambda^{{\bf int}}', [], ...
                ['Ours (N = ', num2str(N), ')'], 'VInt', 'v^{{\bf int}}');   
        end
    end
    
    % ----------------------------------------------------------------------------------------------
    % Plot inner product of interpolated eigenvectors
    % ----------------------------------------------------------------------------------------------
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
PlotRMSE(sPlotParams, vRmseInt, vRmseNys);
