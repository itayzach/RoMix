%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','normal')
%% Illustrate the first eigenfunctions of the 1-D Guassian kernel
% PlotGaussianKernelEigenfunsExample();
%% Load preset
% sPreset = Get1DGridPreset();
% sPreset = Get1DUniformPreset();
% sPreset = Get1DGaussianPreset();
% sPreset = GetSwissRollPreset();
% sPreset = GetBrazilWeatherPreset();
% sPreset = GetTwoMoonsPreset(); %Classifier is okay. Should fix order of VRef relative to V...
sPreset = GetTwoSpiralsPreset();
% sPreset = GetMnistLatentVAEPreset();
%% Parse sPreset
clusterMethod = 'SC'; % 'GMM' / 'SC'
dim                = sPreset.dim;
n                  = sPreset.n;
N                  = sPreset.N;
k                  = sPreset.k;
nnValue            = sPreset.nnValue;
verticesPDF        = sPreset.verticesPDF;
adjacencyType      = sPreset.adjacencyType;
matrixForEigs      = sPreset.matrixForEigs;
nGenDataCompnts    = sPreset.nGenDataCompnts;
dataGenTechnique   = sPreset.dataGenTechnique;
M                  = sPreset.M;
MTilde             = sPreset.MTilde;
omega              = sPreset.omega; % for W
omegaNys           = sPreset.omega; % for Nystrom
omegaRep           = sPreset.omega; % for representer theorem
omegaTilde         = sPreset.omegaTilde; % for our method
gamma1             = sPreset.gamma1;
gamma2             = sPreset.gamma2;
gamma1Rep          = sPreset.gamma1Rep;
gamma2Rep          = sPreset.gamma2Rep;
gmmRegVal          = sPreset.gmmRegVal;
gmmMaxIter         = sPreset.gmmMaxIter;
gmmNumComponents   = sPreset.gmmNumComponents;
sDatasetParams     = sPreset.sDatasetParams;
sDistanceParams    = sPreset.sDistanceParams;
R                  = sPreset.R;
b_debugUseAnalytic = sPreset.b_debugUseAnalytic;
b_forceCtoIdentity = sPreset.b_forceCtoIdentity;
b_normalizePhi     = sPreset.b_normalizePhi;
b_takeEigsFromWRef = sPreset.b_takeEigsFromWRef;
b_flipSign         = sPreset.b_flipSign;
b_pairwiseFlipSign = sPreset.b_pairwiseFlipSign;
b_interpEigenvecs  = sPreset.b_interpEigenvecs;
b_runGraphSignals  = sPreset.b_runGraphSignals;
interpRatio        = N/n;
%% Verify
assert(~b_debugUseAnalytic || (b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
assert(~strcmp(adjacencyType,'NearestNeighbor') || ...
    strcmp(adjacencyType,'NearestNeighbor') && strcmp(verticesPDF,'Grid'))
assert((strcmp(clusterMethod, 'GMM') && n >= dim) || (strcmp(clusterMethod, 'SC'))) 
%% Plot params
sPlotParams = GetPlotParams();
sPlotParams.actualDataDist = verticesPDF;
sPlotParams.matrixForEigs = matrixForEigs;
%% plotInd
if dim == 1
    plotInd = [0,min(4,M-1)];
else
    plotInd = [0,min(8,M-1)];
end
%% Run
mVIntToCompare = zeros(R, N, M);
mVNysToCompare = zeros(R, N, M);
mVRepToCompare = zeros(R, N, M);
mVRefToCompare = zeros(R, N, M);
for r = 1:R
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    fprintf('Iteration r = %d of R = %d\n',r,R)
    sDataset = GenerateDataset(verticesPDF, dim, nGenDataCompnts, n, N, dataGenTechnique, sDatasetParams);
    xTrain   = sDataset.sData.x;
    yTrain   = sDataset.sData.y;
    xInt     = sDataset.sData.xt;
    assert(n == size(xTrain,1));
    assert(N == size(xInt,1));
    if strcmp(sPreset.dataGenTechnique, 'NewPoints')
        warning('sPreset.dataGenTechnique is NewPoints, should make sure its okay...')
        pause(0.5)
    else
        assert(isequal(xTrain, xInt(1:n,:)));
    end

    if r == 1 && sPlotParams.b_plotData && dim <= 3
        PlotDataset(sPlotParams, xTrain, yTrain, 'Training set');
    end
    % ----------------------------------------------------------------------------------------------
    % Original graph
    % ----------------------------------------------------------------------------------------------
    [W, dist, D, Ln] = SimpleCalcAdjacency(xTrain, adjacencyType, sDistanceParams, omega, k, nnValue);    
    [WRef, distRef, DRef, LnRef] = SimpleCalcAdjacency(xInt, adjacencyType, sDistanceParams, omega, k, nnValue);

    if r == 1 && sPlotParams.b_plotWeights
        PlotWeightsMatrix([], W, dist, D, Ln, xTrain, adjacencyType, omega, k);
    end
    % ----------------------------------------------------------------------------------------------
    % Perform numeric eigs
    % ----------------------------------------------------------------------------------------------
    if b_debugUseAnalytic
        warning('review the following {1}...')
        [VRef, adjLambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, ...
            sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
        adjLambdaRef = N*adjLambdaRef;
        matLambda = adjLambdaRef;
    else
        if b_takeEigsFromWRef
            [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, DRef, M, matrixForEigs);

            W = WRef(1:n,1:n);
            dist = distRef(1:n,1:n);
            D = DRef(1:n,1:n);
            V = sqrt(interpRatio)*VRef(1:n,:);
            adjLambda = adjLambdaRef/interpRatio;
            matLambda = matLambdaRef/interpRatio;
            VRefToCompare = VRef;
        else
            [V, adjLambda, matLambda] = EigsByType(W, D, M, matrixForEigs);
        end
    end
    assert(isreal(V), 'V should be real...')
    assert(isreal(matLambda), 'matLambda should be real...')

    % ----------------------------------------------------------------------------------------------
    % Plot numeric V
    % ----------------------------------------------------------------------------------------------
    if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
        figTitle = ['Eigenvectors of ', matrixForEigs];
        if strcmp(adjacencyType, 'GaussianKernel')
            figTitle2 = [' (Gaussian kernel, $\omega = ', num2str(omega), '$'];
        elseif strcmp(adjacencyType, 'NearestNeighbor')
            figTitle2 = [' (k-NN, $k = ', num2str(k), '$'];
        end
        figTitle3 = [', $n = ', num2str(n), '$)'];
        figName = 'V';
        PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
            V, matLambda, '\lambda^{v}', [], [figTitle, figTitle2, figTitle3], figName, 'v');
        figTitle = 'Numeric eigenvalues of ${\bf V}$';
        PlotSpectrum([], [], matLambda, [], [], '\tilde{\lambda}^{v}_m', [], [], figTitle);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Estimate distribution parameters
    % ----------------------------------------------------------------------------------------------
    xTildeTrain = xTrain;
    if strcmp(clusterMethod, 'GMM')
        sDistParams = EstimateDistributionParameters(xTildeTrain, gmmNumComponents, gmmRegVal, gmmMaxIter);
    elseif strcmp(clusterMethod, 'SC')
        sDistParams = EstDistParamsSpectClust(xTildeTrain, W, gmmNumComponents);
    end
    % ----------------------------------------------------------------------------------------------
    % Plot distribution estimation
    % ----------------------------------------------------------------------------------------------
    if r == 1 && sPlotParams.b_plotDataVsGmm && dim <= 3
        nGmmPoints = n;
        pltTitle = ['Dataset with n = ', num2str(n), ' points'];
        if strcmp(clusterMethod, 'GMM')
            plt2Title = ['Generated ' num2str(nGmmPoints), ' points from ' clusterMethod ' with nEstComp = ' num2str(gmmNumComponents)];
        elseif strcmp(clusterMethod, 'SC')
            plt2Title = ['Painted ' num2str(n), ' according to ' clusterMethod ' with nEstComp = ' num2str(gmmNumComponents)];
        end
        windowStyle = 'normal';
        PlotDataset(sPlotParams, xTrain, yTrain, pltTitle, sDistParams.GMModel, nGmmPoints, plt2Title, windowStyle);

    end
    % ----------------------------------------------------------------------------------------------
    % Calculate lambdaAnalyticTilde and PhiTilde(xTildeTrain)
    % ----------------------------------------------------------------------------------------------
    sKernelParams = CalcKernelParams(sDistParams, omegaTilde);
    [sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
         = CalcAnalyticEigenvalues(MTilde, sKernelParams);
    [ PhiTilde, lambdaAnalyticTilde ] = ...
        CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeTrain, b_normalizePhi);
    invLambda = diag(1./lambdaAnalyticTilde);

    % ----------------------------------------------------------------------------------------------
    % Plot PhiTilde
    % ----------------------------------------------------------------------------------------------
    if r == 1 && sPlotParams.b_plotTildeFiguresForDebug && dim <= 3
        figTitle = 'Eigenfunctions of the Gaussian kernel (on $n$ nodes)';
        figName = 'PhiTilde';
        PlotEigenfuncvecScatter([], 'Gaussian', xTildeTrain, [], 0, 4, ...
            PhiTilde, [], [], [], figTitle, figName, '\phi' );
        figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        PlotSpectrum([], [], lambdaAnalyticTilde, [], [], '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
    end
    
    % ----------------------------------------------------------------------------------------------
    % Calculate eigenfunctions values at xInt
    % ----------------------------------------------------------------------------------------------
    [PhiTildeInt, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xInt, b_normalizePhi);

    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = CalcDistance(xTrain, xInt, sDistanceParams);
    WTrainInt = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    lambdaNys = adjLambda*sqrt(interpRatio);
    VNys = WTrainInt.'*V*diag(1./lambdaNys);
    
    if b_interpEigenvecs
        % ----------------------------------------------------------------------------------------------
        % Learn eigenvectors to eigenfunctions transformation (C)
        % ----------------------------------------------------------------------------------------------
        if b_forceCtoIdentity
            C = zeros(MTilde, M);
            C(1:M,1:M) = eye(M);
        else
            b_maskDataTermCMatrix = false;
            C = EigsRLS(PhiTilde, gamma1, gamma2, invLambda, Ln, V, b_maskDataTermCMatrix);
        end
        % ----------------------------------------------------------------------------------------------
        % Interpolate with our method
        % ----------------------------------------------------------------------------------------------
        VInt = (1/sqrt(interpRatio))*PhiTildeInt*C;
        VIntToCompare = VInt;
        % ----------------------------------------------------------------------------------------------
        % Interpolate with Nystrom
        % ----------------------------------------------------------------------------------------------
        if b_flipSign
            VNys = FlipSign(VInt, VNys, b_pairwiseFlipSign);
        end
        VNysToCompare = VNys;
        % ----------------------------------------------------------------------------------------------
        % Learn alpha using the Representer theorem
        % ----------------------------------------------------------------------------------------------
        b_normalizeAlpha = true;
        mAlpha = LapRLS(W, V, Ln, gamma1Rep, gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
        
        if r == 1 && sPlotParams.b_plotC
            CRep = diag(lambdaAnalyticTilde)*PhiTilde.'*mAlpha;
            PlotCoeffsMatrix(C, '${\bf C}$', CRep, ...
                '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
            %             PlotCoefficients(sPlotParams, C(:,1:5), lambdaAnalyticTilde)
        end
        % ----------------------------------------------------------------------------------------------
        % Interpolate with Representer theorem
        % ----------------------------------------------------------------------------------------------
        VRep = WTrainInt.'*mAlpha;
        if b_interpEigenvecs
            if b_flipSign
                VRep = FlipSign(VInt, VRep, b_pairwiseFlipSign);
            end
            VRepToCompare = VRep;
        end
        % ----------------------------------------------------------------------------------------------
        % Calculate reference
        % ----------------------------------------------------------------------------------------------
        if b_debugUseAnalytic
            [PhiRef, adjLambdaRef] = SimpleCalcAnalyticEigenfunctions(xInt, omega, ...
                sDatasetParams.sigma{1}, sDatasetParams.mu{1}, M);
            VRef = PhiRef;
            if b_flipSign
                VRefToCompare = FlipSign(VInt, VRef, b_pairwiseFlipSign);
            end
        elseif ~b_takeEigsFromWRef
            [VRef, adjLambdaRef, matLambdaRef] = EigsByType(WRef, DRef, M, matrixForEigs);
            LnRef = DRef - WRef;
            if b_flipSign
                VRef = FlipSign(VInt, VRef, b_pairwiseFlipSign);
            end
            VRefToCompare = VRef;
        end

        if r == 1 && sPlotParams.b_plotC
            Cint = EigsRLS(PhiTildeInt, gamma1, gamma2, invLambda, LnRef, VRefToCompare, b_maskDataTermCMatrix);
            PlotCoeffsMatrix(C, '${\bf C}$', Cint, '${\bf C^{\bf int}}$');
        end

        if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
            if dim == 1
                PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                    VRefToCompare, [], [], [], ['Reference vs. Representer theorem (N = ', num2str(N), ')'], ...
                    'VRep', 'v^{{\bf ref}}', VRepToCompare, 'v^{{\bf rep}}');
                PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                    VRefToCompare, [], [], [], ['Reference vs. Nystrom (N = ', num2str(N), ')'], ...
                    'VNys', 'v^{{\bf ref}}', VNysToCompare, 'v^{{\bf nys}}');
                PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                    VRefToCompare, [], [], [], ['Reference vs. Ours (N = ', num2str(N), ')'], ...
                    'VInt', 'v^{{\bf ref}}', VIntToCompare, 'v^{{\bf int}}');
            else
                cmap = PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                    VRefToCompare, [], [], [], ...
                    ['Reference (N = ', num2str(N), ')'], 'VRef', 'v^{{\bf ref}}');
                PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, cmap, plotInd(1), plotInd(2), ...
                    VRepToCompare, [], [], [], ...
                    ['Representer theorem (N = ', num2str(N), ')'], 'VRep', 'v^{{\bf rep}}');
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
            PlotInnerProductMatrix([], VRef, [], pltTitle, figName);
            pltTitle = 'VInt - ${\bf V}_{{\bf int}}^T {\bf V}_{{\bf int}}$';
            figName = 'Vint';
            PlotInnerProductMatrix([], VInt, [], pltTitle, figName);
            pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
            figName = 'VNys';
            PlotInnerProductMatrix([], VNys, [], pltTitle, figName);
            pltTitle = 'VRep - ${\bf V}_{{\bf rep}}^T {\bf V}_{{\bf rep}}$';
            figName = 'VRep';
            PlotInnerProductMatrix([], VRep, [], pltTitle, figName);
    %         pltTitle = '$\int \phi_i(x) \phi_j(x) p(x) dx = \Phi^T$diag(Pr)$\Phi$';
    %         figName = 'PhiTilde';
    %         vPrTilde = sDistParams.vPr;
    %         PlotInnerProductMatrix([], PhiTilde, vPrTilde, pltTitle, figName);

        end
    end
    % ----------------------------------------------------------------------------------------------
    % Graph signal
    % ----------------------------------------------------------------------------------------------
    if b_runGraphSignals
        if isfield(sDataset.sData, 'y') && ~isempty(sDataset.sData.y)
            assert(~isempty(sDataset.sData.yt));
            mSig = sDataset.sData.y;
            mSigRef = sDataset.sData.yt;
        else
            [mSig, mSigRef] = GenerateSyntheticGraphSignal(V, VRef);
        end
        % ------------------------------------------------------------------------------------------
        % Coeffs
        % ------------------------------------------------------------------------------------------
        mSigHatV = V'*mSig; % same as pinv(V)*sig...
        if ismember(sPreset.verticesPDF, {'USPS'})  
            mSigMasked = sDataset.sData.ymasked;
            mSigHatPhi = EigsRLS(PhiTilde, gamma1, gamma2, invLambda, Ln, mSigMasked, sPreset.b_maskDataFitTerm);
        else
            mSigHatPhi = EigsRLS(PhiTilde, gamma1, gamma2, invLambda, Ln, mSig, sPreset.b_maskDataFitTerm);
        end
        b_normalizeAlpha = true;
        mAlpha = LapRLS(W, mSig, Ln, gamma1Rep, gamma2Rep, interpRatio, b_normalizeAlpha, sPreset.b_maskDataFitTerm);
        
        % Just for reference
        mSigRefHatPhi = EigsRLS(PhiTildeInt, gamma1, gamma2, invLambda, LnRef, mSigRef, sPreset.b_maskDataFitTerm);
        % ------------------------------------------------------------------------------------------
        % Signals
        % ------------------------------------------------------------------------------------------
        mSigRecV = V*mSigHatV;
        mSigRecPhi = PhiTilde*mSigHatPhi;
        mSigRecRep = W.'*mAlpha;
        mSigInt = PhiTildeInt*mSigHatPhi;
        mSigNys = VNys*mSigHatV;
        mSigRep = WTrainInt.'*mAlpha;

        mSigIntRef = PhiTildeInt*mSigRefHatPhi;
        if (ismember(sPreset.verticesPDF, {'TwoMoons', 'TwoSpirals'})) || ...
                (ismember(sPreset.verticesPDF, {'MnistLatentVAE'}) && isequal(sDatasetParams.vPossibleLabels, [0, 1]))
            mSigRecV = sign(mSigRecV);
            mSigRecPhi = sign(mSigRecPhi);
            mSigRecRep = sign(mSigRecRep);
            mSigInt = sign(mSigInt);
            mSigNys = sign(mSigNys);
            mSigRep = sign(mSigRep);

            mSigIntRef = sign(mSigIntRef);
        elseif ismember(sPreset.verticesPDF, {'MnistLatentVAE'})    
            error('I dont think round is the way to go here...')
            mSigRecV = round(mSigRecV);
            mSigRecPhi = round(mSigRecPhi);
            mSigRecRep = round(mSigRecRep);
            mSigInt = round(mSigInt);
            mSigNys = round(mSigNys);
            mSigRep = round(mSigRep);

            mSigIntRef = round(mSigIntRef);
        elseif ismember(sPreset.verticesPDF, {'USPS', 'MNIST'})  
            [~, mSig] = max(mSig,[],2);
            [~, mSigRef] = max(mSigRef,[],2);

            [~, mSigRecV] = max(mSigRecV,[],2);
            [~, mSigRecPhi] = max(mSigRecPhi,[],2);
            [~, mSigRecRep] = max(mSigRecRep,[],2);
            [~, mSigInt] = max(mSigInt,[],2);
            [~, mSigNys] = max(mSigNys,[],2);
            [~, mSigRep] = max(mSigRep,[],2);

            [~, mSigIntRef] = max(mSigIntRef,[],2);
        end
        
               
        % ------------------------------------------------------------------------------------------
        % Plots
        % ------------------------------------------------------------------------------------------
        sigIndToPlot  = min(7, size(mSig,2));
        vSigHatV      = mSigHatV(:,sigIndToPlot);
        vSigRefHatPhi = mSigRefHatPhi(:,sigIndToPlot);
        vSigHatPhi    = mSigHatPhi(:,sigIndToPlot);
        vAlpha        = mAlpha(:,sigIndToPlot);
        
        vSig       = mSig(:,sigIndToPlot);
        vSigRef    = mSigRef(:,sigIndToPlot);
        vSigRecV   = mSigRecV(:,sigIndToPlot);
        vSigRecPhi = mSigRecPhi(:,sigIndToPlot);
        vSigRecRep = mSigRecRep(:,sigIndToPlot);
        vSigInt    = mSigInt(:,sigIndToPlot);
        vSigIntRef = mSigIntRef(:,sigIndToPlot);
        vSigNys    = mSigNys(:,sigIndToPlot);
        vSigRep    = mSigRep(:,sigIndToPlot);

        if r == 1 && sPlotParams.b_plotC
            CRep = diag(lambdaAnalyticTilde)*PhiTilde.'*mAlpha;
            PlotCoeffsMatrix(mSigHatPhi, '${\bf C}$', CRep, ...
                '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
            PlotCoefficients([], mSigHatPhi(:,sigIndToPlot), lambdaAnalyticTilde)
        end

        
%         PlotGraphSignalAnalysis(vSig, vSigRecPhi, vSigRecV, vSigRef, vSigInt, vSigNys, ...
%             vSigHatPhi, vSigHatV, vSigRefHatPhi);
       
%         PlotGraphSignalErrors(sPlotParams, [1, n], vSig, {vSigRecPhi, vSigRecV, vSigRecRep}, ...
%             {'|s-s_{\Phi}^{{\bf rec}}|', '|s-s_{V}^{{\bf rec}}|', '|s-s_{{\bf K}}^{{\bf rec}}|'},  ...
%             'Projection error ($n$ given nodes)')
%         if n+1 < N
%             PlotGraphSignalErrors(sPlotParams, [n+1, N], vSigRef(n+1:N), ...
%                 {vSigInt(n+1:N), vSigNys(n+1:N), vSigRep(n+1:N)}, ...
%                 {'|s^{{\bf ref}}-s^{{\bf int}}|', '|s^{{\bf ref}}-s^{{\bf nys}}|', '|s^{{\bf ref}}-s^{{\bf rep}}|'},  ...
%                 'Interpolation error ($N-n$ nodes)')
%         end
%         
%         PlotGraphSignalErrors(sPlotParams, [1, N], vSigRef, {vSigInt, vSigNys, vSigRep}, ...
%             {'|s^{{\bf ref}}-s^{{\bf int}}|', '|s^{{\bf ref}}-s^{{\bf nys}}|', '|s^{{\bf ref}}-s^{{\bf rep}}|'},  ...
%             'Total error ($N$ nodes)')
        
        if r == 1 && sPlotParams.b_plotC
            PlotCoefficients(sPlotParams, vSigHatPhi, lambdaAnalyticTilde)        
            PlotCoefficients(sPlotParams, vSigRefHatPhi, lambdaAnalyticTilde)        
        end
        
        if r == 1 && dim <=3 &&  ~(ismember(sPreset.verticesPDF, {'TwoMoons'}))
            PlotGraphSignals(sPlotParams, ['Graph signals on given $n=' num2str(n) '$ nodes (Train set)'], 'TrainSet', ...
                {xTrain, xTrain, xTrain, xTrain}, ...
                {vSig, vSigRecPhi, vSigRecV, vSigRecRep}, ...
                {'$s$', '$s_{\Phi}^{{\bf rec}}$', '$s_{V}^{{\bf rec}}$', '$s_{K}^{{\bf rec}}$'}, ...
                {n, n, n, n});
            cmap = PlotGraphSignals(sPlotParams, ['Graph signals on all $N=' num2str(N) '$ nodes (Train \& Test sets)'], 'TrainAndTestSet', ...
                {xInt, xInt, xInt, xInt}, ...
                {vSigRef, vSigInt, vSigNys, vSigRep}, ...
                {'$s^{{\bf ref}}$', '$s^{{\bf int}}$', '$s^{{\bf nys}}$', '$s^{{\bf rep}}$'}, ...
                {n, n, n, n});
            
            if sPlotParams.b_plotGmmSignal
                % GMM
                nGmmPoints = 1000;
                [xGmm,compIdx] = random(sDistParams.GMModel, nGmmPoints);
                [PhiTildeGmm, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xGmm, b_normalizePhi);
                mSigGmm = PhiTildeGmm*mSigHatPhi;
                vSigGmm = mSigGmm(:,sigIndToPlot);
                if dim == 2
                    xylim(1) = min(xTrain(:,1));
                    xylim(2) = max(xTrain(:,1));
                    xylim(3) = min(xTrain(:,2));
                    xylim(4) = max(xTrain(:,2));
                else
                    xylim = [];
                end
                PlotGraphSignals(sPlotParams, 'GMM Graph signal', 'GMM', {xGmm}, {vSigGmm}, ...
                    {'$s^{{\bf gmm}}$'}, {nGmmPoints}, xylim, cmap);
            end
        end
        
        % ------------------------------------------------------------------------------------------
        % Accuracy
        % ------------------------------------------------------------------------------------------
        if (ismember(sPreset.verticesPDF, {'TwoMoons', 'TwoSpirals'})) || ...
                (ismember(sPreset.verticesPDF, {'MnistLatentVAE'}) && isequal(sDatasetParams.vPossibleLabels, [0, 1]))
            vPredictionsInt = mSigInt;
            vTestIdx = find(mSigRef ~= 0);
            testAccInt = 100*sum(vPredictionsInt(vTestIdx) == mSigRef(vTestIdx)) / length(vTestIdx);
            fprintf('Ours Test accuracy = %.2f%%\n', testAccInt);
            
            sClassifier.MTilde = MTilde;
            sClassifier.omega = omegaTilde;
            sClassifier.gamma_A = gamma1;
            sClassifier.gamma_I = gamma2;
            sClassifier.error = 100 - testAccInt;
            sClassifier.vPhi_xTrain_c = mSigRecPhi;
            xMax = max([xTrain(:); xInt(:)]);
            xMin = min([xTrain(:); xInt(:)]);
            step = (xMax - xMin)/100;
            x1 = xMin:step:xMax;
            x2 = x1;
            [sClassifier.XX1,sClassifier.XX2] = meshgrid(x1,x2);
            xMeshGrid = [sClassifier.XX1(:) sClassifier.XX2(:)];   
            PhiTildeMeshgrid = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xMeshGrid, b_normalizePhi);
            sClassifier.mPhi_X_c = reshape(PhiTildeMeshgrid*mSigHatPhi, length(x1), length(x2));
            PlotGraphSignalAnalysis(mSig, mSigRecPhi, mSigRecV, mSigRef, mSigInt, mSigNys, ...
                mSigHatPhi, mSigHatV, mSigRefHatPhi);
        
            PlotClassifier(sPlotParams, sDataset, sClassifier);


            vPredictionsRep = mSigRep;
            testAccRepThm = 100*sum(vPredictionsRep(vTestIdx) == mSigRef(vTestIdx)) / length(vTestIdx);
            fprintf('Rep. Thm Test accuracy = %.2f%%\n', testAccRepThm);

            sLapRLSClassifier.alpha = mAlpha;
            sLapRLSClassifier.gammas = [gamma1Rep, gamma2Rep];
            sLapRLSClassifier.test_error = 100 - testAccRepThm;
            PlotLapRLS(sPlotParams, sDataset, sDistanceParams, omega, sLapRLSClassifier);

        elseif ismember(sPreset.verticesPDF, {'MnistLatentVAE', 'USPS', 'MNIST'})
            vPredictionsRecRep = mSigRecRep;
            testAccRecRep = 100*sum(vPredictionsRecRep == mSig) / n;
            fprintf('\nRep Thm train accuracy = %.2f%%\n', testAccRecRep);
            vPredictionsRep = mSigRep;
            testAccRep = 100*sum(vPredictionsRep == mSigRef) / N;
            fprintf('Rep Thm test accuracy  = %.2f%%\n', testAccRep);
            
            vPredictionsRecV = mSigRecV;
            testAccRecV = 100*sum(vPredictionsRecV == mSig) / n;
            fprintf('\nNystrom train accuracy = %.2f%%\n', testAccRecV);
            vPredictionsNys = mSigNys;
            testAccNys = 100*sum(vPredictionsNys == mSigRef) / N;
            fprintf('Nystrom test accuracy  = %.2f%%\n', testAccNys);
            
            vPredictionsRecPhi = mSigRecPhi;
            testAccRecPhi = 100*sum(vPredictionsRecPhi == mSig) / n;
            fprintf('\nOurs train accuracy    = %.2f%%\n', testAccRecPhi);
            vPredictionsInt = mSigInt;
            testAccInt = 100*sum(vPredictionsInt == mSigRef) / N;
            fprintf('Ours test accuracy     = %.2f%%\n', testAccInt);
            
            vPredictionsIntRef = mSigIntRef;
            testAccIntRef = 100*sum(vPredictionsIntRef == mSigRef) / N;
            fprintf('Ours test accuracy     = %.2f%%\n', testAccIntRef);
        else
            [vRmseRecRep, vMseRecRep, vAccRecRep, mErrRecRep, vCohRecRep] = ...
                CalcErrAndAcc(mSigRecRep, mSig, 'Representer (train)');
            [vRmseRep, vMseRep, vAccRep, mErrRep, vCohRep] = ...
                CalcErrAndAcc(mSigRep, mSigRef, 'Representer (test)');
            [vRmseRecV, vMseRecV, vAccRecV, mErrRecV, vCohRecV] = ...
                CalcErrAndAcc(mSigRecV, mSig, 'Nystrom (train)');
            [vRmseNys, vMseNys, vAccNys, mErrNys, vCohNys] = ...
                CalcErrAndAcc(mSigNys, mSigRef, 'Nystrom (test)');
            [vRmseRecPhi, vMseRecPhi, vAccRecPhi, mErrRecPhi, vCohRecPhi] = ...
                CalcErrAndAcc(mSigRecPhi, mSig, 'EigsRLS (train)');
            [vRmseInt, vMseInt, vAccInt, mErrInt, vCohInt] = ...
                CalcErrAndAcc(mSigInt, mSigRef, 'EigsRLS (test)');
            
%             PlotRMSE(sPlotParams, [vRmseRecPhi, vRmseRecV, vRmseRecRep], ...
%                 {'RMSE$(f^{{\bf int}}_m, f^{{\bf ref}}_m)$', 'RMSE$(f^{{\bf nys}}_m, f^{{\bf ref}}_m)$', ...
%                 'RMSE$(f^{{\bf rep}}_m, f^{{\bf ref}}_m)$'}, 'Train RMSE');
            PlotAccuracy(sPlotParams, [vAccRecPhi, vAccRecV, vAccRecRep], ...
                {'Acc$(f^{{\bf int}}_m, f^{{\bf ref}}_m)$', 'Acc$(f^{{\bf nys}}_m, f^{{\bf ref}}_m)$', ...
                'Acc$(f^{{\bf rep}}_m, f^{{\bf ref}}_m)$'}, 'TrainAcc', [], sDatasetParams.monthNames, 'Train accuracy');
            PlotEvalMetric(sPlotParams, [vCohRecPhi, vCohRecV, vCohRecRep], ...
                {'Coh$(s_{\Phi}^{{\bf rec}}, s)$', 'Coh$(s_{V}^{{\bf rec}}, s)$', ...
                'Coh$(s_{K}^{{\bf rec}}, s)$'}, 'Coh_graph_signals_', [], sDatasetParams.monthNames, 'Train coherence');
%             PlotEvalMetric(sPlotParams, [vMseRecPhi, vMseRecV, vMseRecRep], ...
%                 {'MSE$(s_{\Phi}^{{\bf rec}}, s)$', 'MSE$(s_{V}^{{\bf rec}}, s)$', ...
%                 'MSE$(s_{K}^{{\bf rec}}, s)$'}, 'Mse_graph_signals_', [], sDatasetParams.monthNames, 'Train MSE');

%             PlotRMSE(sPlotParams, [vRmseInt, vRmseNys, vRmseRep], ...
%                 {'RMSE$(f^{{\bf int}}_m, f^{{\bf ref}}_m)$', 'RMSE$(f^{{\bf nys}}_m, f^{{\bf ref}}_m)$', ...
%                 'RMSE$(f^{{\bf rep}}_m, f^{{\bf ref}}_m)$'}, 'Test RMSE');
            PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], ...
                {'Acc$(f^{{\bf int}}_m, f^{{\bf ref}}_m)$', 'Acc$(f^{{\bf nys}}_m, f^{{\bf ref}}_m)$', ...
                'Acc$(f^{{\bf rep}}_m, f^{{\bf ref}}_m)$'}, 'TestAcc', [], sDatasetParams.monthNames, 'Train \& Test accuracy');
            PlotEvalMetric(sPlotParams, [vCohInt, vCohNys, vCohRep], ...
                {'Coh$(f^{{\bf int}}_m, f^{{\bf ref}}_m)$', 'Coh$(f^{{\bf nys}}_m, f^{{\bf ref}}_m)$', ...
                'Coh$(f^{{\bf rep}}_m, f^{{\bf ref}}_m)$'}, 'Coh_graph_signals_', [], sDatasetParams.monthNames, 'Train \& Test coherence');
%             PlotEvalMetric(sPlotParams, [vMseInt, vMseNys, vMseRep], ...
%                 {'MSE$(f^{{\bf int}}_m, f^{{\bf ref}}_m)$', 'MSE$(f^{{\bf nys}}_m, f^{{\bf ref}}_m)$', ...
%                 'MSE$(f^{{\bf rep}}_m, f^{{\bf ref}}_m)$'}, 'Mse_graph_signals_', [], sDatasetParams.monthNames, 'Train \& Test  MSE');
        end
        
    end
    
    % Save
    if b_interpEigenvecs
        mVIntToCompare(r,:,:) = VIntToCompare;
        mVNysToCompare(r,:,:) = VNysToCompare;
        mVRepToCompare(r,:,:) = VRepToCompare;
        mVRefToCompare(r,:,:) = VRefToCompare;
    end
    
end
if b_interpEigenvecs
    [vRmseInt, vMseInt, vAccInt, mErrInt, vCohInt] = CalcErrAndAcc(mVIntToCompare, mVRefToCompare, 'Analytic');
    [vRmseNys, vMseNys, vAccNys, mErrNys, vCohNys] = CalcErrAndAcc(mVNysToCompare, mVRefToCompare, 'Nystrom');
    [vRmseRep, vMseRep, vAccRep, mErrRep, vCohRep] = CalcErrAndAcc(mVRepToCompare, mVRefToCompare, 'Representer');
%     PlotRMSE(sPlotParams, [vRmseInt, vRmseNys, vRmseRep], ...
%         {'RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
%          'RMSE$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'});
    PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], ...
        {'Acc$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'Acc$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
         'Acc$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'}, ['Acc_eigs_0_to_' num2str(M-1)]);
    PlotEvalMetric(sPlotParams, [vCohInt, vCohNys, vCohRep], ...
        {'Coh$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'Coh$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
         'Coh$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'}, ['Coh_eigs_0_to_' num2str(M-1)]);
%     PlotEvalMetric(sPlotParams, [vMseInt, vMseNys, vMseRep], ...
%         {'MSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'MSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
%          'MSE$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'}, ['Mse_eigs_0_to_' num2str(M-1)]);
end
