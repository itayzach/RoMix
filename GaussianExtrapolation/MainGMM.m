%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','docked')
%% Plot params
sPlotParams = GetPlotParams();
%% Dataset parameters
dim                 = 1;
n                   = 249;
N                   = 299;
k                   = round(0.01*N);
nnValue             = 'ZeroOne'; % 'ZeroOne' / 'Distance'
verticesPDF         = 'BrazilWeather'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'SwissRoll' / 'MnistLatentVAE' / 'CoraLatentVGAE' / 'BrazilWeather'
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
MTilde             = 50;
%% GMM params
gmmRegVal          = 1e-3;
gmmMaxIter         = 2000;
gmmNumComponents   = 3;
%% Method parameters
b_debugUseAnalytic = false;
b_forceCtoIdentity = false;
b_normalizePhi     = false;
b_takeEigsFromWRef = false;
%% Verify
assert(~b_debugUseAnalytic || (b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
assert(~strcmp(adjacencyType,'NearestNeighbor') || ...
    strcmp(adjacencyType,'NearestNeighbor') && strcmp(verticesPDF,'Grid'))
if strcmp(verticesPDF, 'BrazilWeather')
    distType = 'Haversine'; % 'Euclidean' / 'Haversine'
else
    distType = 'Euclidean'; % 'Euclidean' / 'Haversine'
end

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
assert(isequal(sDataset.sData.x, sDataset.sData.xt(1:n,:)));
%% plotInd
if dim == 1
    plotInd = [0,min(4,M-1)];
else
    plotInd = [0,min(11,M-1)];
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
        [VRef, adjLambdaRef] = SimpleCalcAnalyticEigenfunctions(xTrainInt, omega, ...
            sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        adjLambdaRef = N*adjLambdaRef;
        matLambda = adjLambdaRef;
    else
        if b_takeEigsFromWRef
            [WRef, distRef, DRef] = SimpleCalcAdjacency(xInt, adjacencyType, distType, omega, k, nnValue);
            [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, DRef, M, matrixForEigs);

            W = WRef(1:n,1:n);
            dist = distRef(1:n,1:n);
            D = DRef(1:n,1:n);
            V = VRef(1:n,:);
            adjLambda = adjLambdaRef;
            matLambda = matLambdaRef;
            VRefToCompare = VRef;
        else
            [W, dist, D] = SimpleCalcAdjacency(xTrain, adjacencyType, distType, omega, k, nnValue);
            [V, adjLambda, matLambda] = EigsByType(W, D, M, matrixForEigs);
            V = V/sqrt(interpRatio);
        end
        
        if r == 1 && sPlotParams.b_plotWeights
            PlotWeightsMatrix(sPlotParams, W, dist, D, xTrain, adjacencyType, verticesPDF, omega, k);
        end
        assert(isreal(V), 'V should be real...')
        assert(isreal(matLambda), 'matLambda should be real...')
    end

    % ----------------------------------------------------------------------------------------------
    % Calculate lambdaAnalyticTilde and PhiTilde(xTildeTrain)
    % ----------------------------------------------------------------------------------------------
    xTildeTrain = xTrain;
    sDistParams = EstimateDistributionParameters(xTildeTrain, gmmNumComponents, gmmRegVal, gmmMaxIter);
    vPrTilde = sDistParams.vPr;
    sKernelParams = GetKernelParams(sDistParams, omegaTilde);
    [sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
        = CalcAnalyticEigenvalues(MTilde, sKernelParams, dim, gmmNumComponents);
    [ PhiTilde, lambdaAnalyticTilde ] = ...
        CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeTrain, b_normalizePhi);
    
    if r == 1 && sPlotParams.b_plotDataVsGmm && dim <= 3
        nGmmPoints = n;
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
%         C = n^(dim-1) * (PhiTilde.' * diag(vPrTilde)* V);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Plot PhiTilde
    % ----------------------------------------------------------------------------------------------
    if r == 1 && sPlotParams.b_plotTildeFiguresForDebug && dim <= 3
        figTitle = 'Eigenfunctions of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        figName = 'PhiTilde';
        PlotEigenfuncvecScatter(sPlotParams, 'Gaussian', xTildeTrain, [], 0, 4, ...
            PhiTilde, [], [], [], figTitle, figName, '\tilde{\phi}' );
        figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        PlotSpectrum(sPlotParams, sDataset, [], lambdaAnalyticTilde, [], [], ...
            '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
    end
    PhiTildeCondNum = cond(PhiTilde);
    recommendedMTilde = MTilde;
    while(PhiTildeCondNum > 1e3 && recommendedMTilde >= 10)
        fprintf('Condition number %.3f is too large. Checking with MTilde = %d\n', PhiTildeCondNum, recommendedMTilde)
        recommendedMTilde = recommendedMTilde - 10;
        PhiTildeCondNum = cond(PhiTilde(1:recommendedMTilde,:));
    end
    if MTilde ~= recommendedMTilde
        warning('Consider MTilde = %d with condition number %.3f', recommendedMTilde, PhiTildeCondNum);
        pause(1);
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
        cmap = PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            V, matLambda, '\lambda^{v}', [], [figTitle, figTitle2, figTitle3], figName, 'v');
        figTitle = 'Numeric eigenvalues of ${\bf V}$';
        PlotSpectrum(sPlotParams, sDataset, [], matLambda, [], [], ...
                '\tilde{\lambda}^{v}_m', [], [], figTitle);
     end

    % ----------------------------------------------------------------------------------------------
    % Interpolate with our method
    % ----------------------------------------------------------------------------------------------
    xTildeInt = xInt;
    [PhiTildeInt, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeInt, b_normalizePhi);
    PhiTildeIntCondNum = cond(PhiTildeInt);
%     assert(PhiTildeIntCondNum < 1e3, 'Condition number is too large...');
    VInt = PhiTildeInt*C;
    VIntToCompare = VInt;
       
    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = CalcDistance(xTrain, xInt, distType);
    B = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    VNys = B.'*V*diag(1./adjLambda);
    VNys = FlipSign(VInt, VNys);
    VNysToCompare = VNys;
    lambdaNys = adjLambda;
    
    % ----------------------------------------------------------------------------------------------
    % Calculate reference
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        [PhiRef, adjLambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, ...
            sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        VRef = PhiRef;
        VRefToCompare = FlipSign(VInt, VRef);
    elseif ~b_takeEigsFromWRef
        [WRef, ~, DRef] = SimpleCalcAdjacency(xInt, adjacencyType, distType, omega, k, nnValue);
        [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, DRef, M, matrixForEigs);
        VRef = FlipSign(VInt, VRef);
        VRefToCompare = VRef;
    end

    if r == 1 && sPlotParams.b_plotC
        Cint = pinv(PhiTildeInt)*VRefToCompare;
        PlotCoeffsMatrix(C, Cint);
    end
    
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        if dim == 1
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VRefToCompare, [], [], [], ['Reference vs. Nystrom (N = ', num2str(N), ')'], ...
                'VRef', 'v^{{\bf ref}}', VNysToCompare, 'v^{{\bf nys}}');
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                VRefToCompare, [], [], [], ['Reference vs. Ours (N = ', num2str(N), ')'], ...
                'VRef', 'v^{{\bf ref}}', VIntToCompare, 'v^{{\bf int}}');
        else
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
                VRefToCompare, [], [], [], ...
                ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf ref}}');
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
                VNysToCompare, [], [], [], ...
                ['Nystrom (N = ', num2str(N), ')'], 'VNys', 'v^{{\bf nys}}');
            PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
                VIntToCompare, [], [], [], ...
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
%         pltTitle = '$\int \phi_i(x) \phi_j(x) p(x) dx = \Phi^T$diag(Pr)$\Phi$';
%         figName = 'PhiTilde';
%         PlotInnerProductMatrix(sPlotParams, dim, vPrTilde, 'IP_Matrix', [], PhiTilde, pltTitle, figName);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Graph signal
    % ----------------------------------------------------------------------------------------------
    if r == 1
        if isfield(sDataset, 'graphSignal')
            sig = sDataset.graphSignal;
            sigRef = sDataset.graphSignalInt;
        else
            [sig, sigRef] = GenerateSyntheticGraphSignal(V);
        end

        sigHatV = interpRatio*V'*sig; % same as pinv(V)*sig...
        sigRecV = V*sigHatV;

        sigHatPhi = pinv(PhiTilde)*sig;
        sigRefHatPhi = pinv(PhiTildeInt)*sigRef;
        sigRecPhi = PhiTilde*sigHatPhi;

        sigInt = PhiTildeInt*sigHatPhi;
        sigNys = VNys*sigHatV;

        PlotGraphSignalAnalysis(sig, sigRecPhi, sigRecV, sigRef, sigInt, sigNys, ...
            sigHatPhi, sigHatV, sigRefHatPhi);
        
        PlotGraphSignals(sPlotParams, 'interpolated graph signal', ...
            {xTrain, xTrain, xTrain, xInt, xInt, xInt}, ...
            {sig, sigRecPhi, sigRecV sigRef, sigInt, sigNys}, ...
            {'$s$ ($n$ nodes)', '$s_{\Phi}^{{\bf rec}}$ ($n$ nodes)', '$s_{V}^{{\bf rec}}$ ($n$ nodes)', ...
             '$s^{{\bf ref}}$ ($N$ nodes)', '$s^{{\bf int}}$ ($N$ nodes)', '$s^{{\bf nys}}$ ($N$ nodes)'}, ...
            {n, n, n, n, n, n});
    end
    
    % Save
    mVIntToCompare(r,:,:) = VIntToCompare;
    mVNysToCompare(r,:,:) = VNysToCompare;
    mVRefToCompare(r,:,:) = VRefToCompare;
    
end
vRmseInt = CalcRMSE(mVIntToCompare, mVRefToCompare, 'Analytic');
vRmseNys = CalcRMSE(mVNysToCompare, mVRefToCompare, 'Nystrom');
PlotRMSE(sPlotParams, vRmseInt, vRmseNys);


nGmmPoints = 50000;
[xGmm,compIdx] = random(sDistParams.GMModel, nGmmPoints);
[PhiTildeGmm, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xGmm, b_normalizePhi);
sigGmm = PhiTildeGmm*sigHatPhi;
PlotGraphSignals(sPlotParams, 'interpolated graph signal', ...
    {xGmm}, {sigGmm}, {'$s^{{\bf gmm}}$'}, {n});
