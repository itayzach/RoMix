%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','normal')
%% Load preset
% sPreset = Get1DGridPreset();
% sPreset = GetSwissRollPreset();
sPreset = GetBrazilWeatherPreset();
% sPreset = GetTwoMoonsPreset();
% sPreset = GetTwoSpiralsPreset();
% sPreset = GetMnistLatentVAEPreset();
%% Parse sPreset
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
omega              = sPreset.omega;
omegaNys           = sPreset.omega;
omegaTilde         = sPreset.omegaTilde;
gamma1             = sPreset.gamma1;
gamma2             = sPreset.gamma2;
gmmRegVal          = sPreset.gmmRegVal;
gmmMaxIter         = sPreset.gmmMaxIter;
gmmNumComponents   = sPreset.gmmNumComponents;
b_debugUseAnalytic = sPreset.b_debugUseAnalytic;
b_forceCtoIdentity = sPreset.b_forceCtoIdentity;
b_normalizePhi     = sPreset.b_normalizePhi;
b_takeEigsFromWRef = sPreset.b_takeEigsFromWRef;
b_flipSign         = sPreset.b_flipSign;
b_pairwiseFlipSign = sPreset.b_pairwiseFlipSign;
sDatasetParams     = sPreset.sDatasetParams;
distType           = sPreset.distType;
b_interpEigenvecs  = sPreset.b_interpEigenvecs;
b_runGraphSignals  = sPreset.b_runGraphSignals;
%% Verify
assert(~b_debugUseAnalytic || (b_debugUseAnalytic && strcmp(verticesPDF,'Gaussian')))
assert(~strcmp(adjacencyType,'NearestNeighbor') || ...
    strcmp(adjacencyType,'NearestNeighbor') && strcmp(verticesPDF,'Grid'))
%% Dataset parameters
sDataset = GenerateDataset(verticesPDF, dim, nGenDataCompnts, n, N, dataGenTechnique, sDatasetParams);
rng('default');
dim            = sDataset.dim;
n              = length(sDataset.sData.x);
N              = length(sDataset.sData.xt);
interpRatio    = N/n;
if strcmp(sPreset.dataGenTechnique, 'NewPoints')
    warning('sPreset.dataGenTechnique is NewPoints, should make sure its okay...')
    pause(0.5)
else
    assert(isequal(sDataset.sData.x, sDataset.sData.xt(1:n,:)));
end
%% Plot params
sPlotParams = GetPlotParams();
sPlotParams.sDataset = sDataset;
sPlotParams.matrixForEigs = matrixForEigs;
%% plotInd
if dim == 1
    plotInd = [0,min(4,M-1)];
else
    plotInd = [0,min(8,M-1)];
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
    sDataset = GenerateDataset(verticesPDF, dim, nGenDataCompnts, n, N, dataGenTechnique, sDatasetParams);
    xTrain   = sDataset.sData.x;
    yTrain   = sDataset.sData.y;
    xInt     = sDataset.sData.xt;
    if r == 1 && sPlotParams.b_plotData && dim <= 3
        PlotDataset(sPlotParams, xTrain, yTrain, 'Training set');
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
        [WRef, distRef, DRef, DsqrtInvRef] = SimpleCalcAdjacency(xInt, adjacencyType, distType, omega, k, nnValue);
        LnRef = eye(N) - DsqrtInvRef*WRef*DsqrtInvRef;

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
            [W, dist, D, DsqrtInv] = SimpleCalcAdjacency(xTrain, adjacencyType, distType, omega, k, nnValue);
            [V, adjLambda, matLambda] = EigsByType(W, D, M, matrixForEigs);
            Ln = eye(n) - DsqrtInv*W*DsqrtInv;
%             V = V/sqrt(interpRatio);
        end
        
        if r == 1 && sPlotParams.b_plotWeights
            PlotWeightsMatrix(sPlotParams, W, dist, D, Ln, xTrain, adjacencyType, omega, k);
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
    invLambda = diag(1./lambdaAnalyticTilde);
    if r == 1 && sPlotParams.b_plotDataVsGmm && dim <= 3
        nGmmPoints = n;
        pltTitle = ['Dataset with n = ', num2str(n), ' points'];
        plt2Title = ['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(gmmNumComponents)];

        windowStyle = 'normal';
        PlotDataset(sPlotParams, xTrain, yTrain, pltTitle, sDistParams.GMModel, nGmmPoints, plt2Title, windowStyle);

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
    %         C = n^(dim-1) * (PhiTilde.' * diag(vPrTilde)* V);
        end

        PhiTildeCondNum = cond((PhiTilde).'*(PhiTilde) + gamma1*invLambda);
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
                figTitle2 = [' (Gaussian kernel, $\omega = ', num2str(omega), '$'];
            elseif strcmp(adjacencyType, 'NearestNeighbor')
                figTitle2 = [' (k-NN, $k = ', num2str(k), '$'];
            end
            figTitle3 = [', $n = ', num2str(n), '$)'];
            figName = 'V';
            cmap = PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], plotInd(1), plotInd(end), ...
                V, matLambda, '\lambda^{v}', [], [figTitle, figTitle2, figTitle3], figName, 'v');
            figTitle = 'Numeric eigenvalues of ${\bf V}$';
            PlotSpectrum(sPlotParams, sDataset, [], matLambda, [], [], ...
                    '\tilde{\lambda}^{v}_m', [], [], figTitle);
        end
    end
    % ----------------------------------------------------------------------------------------------
    % Interpolate with our method
    % ----------------------------------------------------------------------------------------------
    xTildeInt = xInt;
    [PhiTildeInt, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTildeInt, b_normalizePhi);
    PhiTildeIntCondNum = cond(PhiTildeInt);
%     assert(PhiTildeIntCondNum < 1e3, 'Condition number is too large...');
    if b_interpEigenvecs
        VInt = (1/sqrt(interpRatio))*PhiTildeInt*C;
        VIntToCompare = VInt;
    end
       
    % ----------------------------------------------------------------------------------------------
    % Interpolate with Nystrom
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = CalcDistance(xTrain, xInt, distType);
    B = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    lambdaNys = adjLambda*sqrt(interpRatio);
    VNys = B.'*V*diag(1./lambdaNys);
    if b_interpEigenvecs
        if b_flipSign
            VNys = FlipSign(VInt, VNys, b_pairwiseFlipSign);
        end
        VNysToCompare = VNys;
    end
    % ----------------------------------------------------------------------------------------------
    % Calculate reference
    % ----------------------------------------------------------------------------------------------
    if b_interpEigenvecs
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
            PlotCoeffsMatrix(C, Cint);
        end

        if r == 1 && (sPlotParams.b_plotOrigVsInterpEvecs || sPlotParams.b_plotAllEvecs) && dim <= 3
            if dim == 1
                PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                    VRefToCompare, [], [], [], ['Reference vs. Nystrom (N = ', num2str(N), ')'], ...
                    'VNys', 'v^{{\bf ref}}', VNysToCompare, 'v^{{\bf nys}}');
                PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xInt, [], plotInd(1), plotInd(2), ...
                    VRefToCompare, [], [], [], ['Reference vs. Ours (N = ', num2str(N), ')'], ...
                    'VInt', 'v^{{\bf ref}}', VIntToCompare, 'v^{{\bf int}}');
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
            PlotInnerProductMatrix(VRef, [], sPlotParams, pltTitle, figName);
            pltTitle = 'VInt - ${\bf V}_{{\bf int}}^T {\bf V}_{{\bf int}}$';
            figName = 'Vint';
            PlotInnerProductMatrix(VInt, [], sPlotParams, pltTitle, figName);
            pltTitle = 'VNys - ${\bf V}_{{\bf nys}}^T {\bf V}_{{\bf nys}}$';
            figName = 'VNys';
            PlotInnerProductMatrix(VNys, [], sPlotParams, pltTitle, figName);
    %         pltTitle = '$\int \phi_i(x) \phi_j(x) p(x) dx = \Phi^T$diag(Pr)$\Phi$';
    %         figName = 'PhiTilde';
    %         PlotInnerProductMatrix(PhiTilde, vPrTilde, sPlotParams, pltTitle, figName);

        end
    end
    % ----------------------------------------------------------------------------------------------
    % Graph signal
    % ----------------------------------------------------------------------------------------------
    if r == 1 && b_runGraphSignals
        if isfield(sDataset.sData, 'y') && ~isempty(sDataset.sData.y)
            assert(~isempty(sDataset.sData.yt));
            sig = sDataset.sData.y;
            sigRef = sDataset.sData.yt;
        else
            [sig, sigRef] = GenerateSyntheticGraphSignal(V, VRef);
        end
      
        % Coeffs
        sigHatV = V'*sig; % same as pinv(V)*sig...
        sigHatPhi = EigsRLS(PhiTilde, gamma1, gamma2, invLambda, Ln, sig, sPreset.b_maskDataFitTerm);
        sigRefHatPhi = EigsRLS(PhiTildeInt, gamma1, gamma2, invLambda, LnRef, sigRef, sPreset.b_maskDataFitTerm);

        % Signals
        sigRecV = V*sigHatV;
        sigRecPhi = PhiTilde*sigHatPhi;
        sigInt = PhiTildeInt*sigHatPhi;
        sigNys = VNys*sigHatV;
        if (ismember(sPreset.verticesPDF, {'TwoMoons', 'TwoSpirals'})) || ...
                (ismember(sPreset.verticesPDF, {'MnistLatentVAE'}) && isequal(sDatasetParams.vPossibleLabels, [0, 1]))
            sigRecV = sign(sigRecV);
            sigRecPhi = sign(sigRecPhi);
            sigInt = sign(sigInt);
            sigNys = sign(sigNys);
        elseif ismember(sPreset.verticesPDF, {'MnistLatentVAE'})    
            sigRecV = round(sigRecV);
            sigRecPhi = round(sigRecPhi);
            sigInt = round(sigInt);
            sigNys = round(sigNys);
        end
        PlotGraphSignalAnalysis(sig, sigRecPhi, sigRecV, sigRef, sigInt, sigNys, ...
            sigHatPhi, sigHatV, sigRefHatPhi);
       
        PlotGraphSignalErrors(sPlotParams, [1, n], sig, sigRecPhi, sigRecV, '|s-s_{\Phi}^{{\bf rec}}|', ...
            '|s-s_{V}^{{\bf rec}}|',  'Projection error ($n$ given nodes)')
        if n+1 < N
            PlotGraphSignalErrors(sPlotParams, [n+1, N], sigRef(n+1:N), sigInt(n+1:N), sigNys(n+1:N), '|s^{{\bf ref}}-s^{{\bf int}}|', ...
                '|s^{{\bf ref}}-s_{{\bf nys}}|',  'Interpolation error ($N-n$ nodes)')
        end
        
        PlotGraphSignalErrors(sPlotParams, [1, N], sigRef, sigInt, sigNys, '|s^{{\bf ref}}-s^{{\bf int}}|', ...
            '|s^{{\bf ref}}-s_{{\bf nys}}|',  'Total error ($N$ nodes)')
        
        cmap = PlotGraphSignals(sPlotParams, 'Graph signals', 'RefIntNys', ...
            {xTrain, xTrain, xTrain, xInt, xInt, xInt}, ...
            {sig, sigRecPhi, sigRecV sigRef, sigInt, sigNys}, ...
            {'$s$ ($n$ nodes)', '$s_{\Phi}^{{\bf rec}}$ ($n$ nodes)', '$s_{V}^{{\bf rec}}$ ($n$ nodes)', ...
             '$s^{{\bf ref}}$ ($N$ nodes)', '$s^{{\bf int}}$ ($N$ nodes)', '$s^{{\bf nys}}$ ($N$ nodes)'}, ...
            {n, n, n, n, n, n});
        
        if (ismember(sPreset.verticesPDF, {'TwoMoons', 'TwoSpirals'})) || ...
                (ismember(sPreset.verticesPDF, {'MnistLatentVAE'}) && isequal(sDatasetParams.vPossibleLabels, [0, 1]))
            vPredictionsInt = sigInt;
            vTestIdx = find(sDataset.sData.yt ~= 0);
            testAccInt = 100*sum(vPredictionsInt(vTestIdx) == sDataset.sData.yt(vTestIdx)) / length(vTestIdx);
            fprintf('Test accuracy = %.2f%%\n', testAccInt);
            
            sClassifier.MTilde = MTilde;
            sClassifier.omega = omegaTilde;
            sClassifier.gamma_A = gamma1;
            sClassifier.gamma_I = gamma2;
            sClassifier.error = 100 - testAccInt;
            sClassifier.vPhi_xTrain_c = sigRecPhi;
            xMax = max([sDataset.sData.x(:); sDataset.sData.xt(:)]);
            xMin = min([sDataset.sData.x(:); sDataset.sData.xt(:)]);
            step = (xMax - xMin)/100;
            x1 = xMin:step:xMax;
            x2 = x1;
            [sClassifier.XX1,sClassifier.XX2] = meshgrid(x1,x2);
            xMeshGrid = [sClassifier.XX1(:) sClassifier.XX2(:)];   
            PhiTildeMeshgrid = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xMeshGrid, b_normalizePhi);
            sClassifier.mPhi_X_c = reshape(PhiTildeMeshgrid*sigHatPhi, length(x1), length(x2));
            PlotGraphSignalAnalysis(sig, sigRecPhi, sigRecV, sigRef, sigInt, sigNys, ...
                sigHatPhi, sigHatV, sigRefHatPhi);
        
            PlotClassifier(sPlotParams, sDataset, sClassifier)
        elseif ismember(sPreset.verticesPDF, {'MnistLatentVAE'})
            vPredictionsRecV = sigRecV;
            testAccRecV = 100*sum(vPredictionsRecV == sDataset.sData.y) / n;
            fprintf('\nNystrom train accuracy = %.2f%%\n', testAccRecV);
            vPredictionsNys = sigNys;
            testAccNys = 100*sum(vPredictionsNys == sDataset.sData.yt) / N;
            fprintf('Nystrom test accuracy  = %.2f%%\n', testAccNys);
            
            vPredictionsRecPhi = sigRecPhi;
            testAccRecPhi = 100*sum(vPredictionsRecPhi == sDataset.sData.y) / n;
            fprintf('Ours train accuracy    = %.2f%%\n', testAccRecPhi);
            vPredictionsInt = sigInt;
            testAccInt = 100*sum(vPredictionsInt == sDataset.sData.yt) / N;
            fprintf('Ours test accuracy     = %.2f%%\n', testAccInt);
        elseif ismember(sPreset.verticesPDF, {'BrazilWeather'})
            vPredictionsRecV = sigRecV;
            testAccRecV = 100*(1 -norm(vPredictionsRecV - sDataset.sData.y) / norm(sDataset.sData.y));
            fprintf('\nNystrom train accuracy = %.2f%%\n', testAccRecV);
            vPredictionsNys = sigNys;
            testAccNys = 100*(1 -norm(vPredictionsNys - sDataset.sData.yt) / norm(sDataset.sData.yt));
            fprintf('Nystrom test accuracy  = %.2f%%\n', testAccNys);
            
            vPredictionsRecPhi = sigRecPhi;
            testAccRecPhi = 100*(1 - norm(vPredictionsRecPhi - sDataset.sData.y) / norm(sDataset.sData.y));
            fprintf('Ours train accuracy    = %.2f%%\n', testAccRecPhi);
            vPredictionsInt = sigInt;
            testAccInt = 100*(1 -norm(vPredictionsInt - sDataset.sData.yt) / norm(sDataset.sData.yt));
            fprintf('Ours test accuracy     = %.2f%%\n', testAccInt);            
        end
        nGmmPoints = 1000;
        [xGmm,compIdx] = random(sDistParams.GMModel, nGmmPoints);
        [PhiTildeGmm, ~] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xGmm, b_normalizePhi);
        sigGmm = PhiTildeGmm*sigHatPhi;
        xylim(1) = min(xTrain(:,1));
        xylim(2) = max(xTrain(:,1));
        xylim(3) = min(xTrain(:,2));
        xylim(4) = max(xTrain(:,2));
        PlotGraphSignals(sPlotParams, 'GMM Graph signal', 'GMM', {xGmm}, {sigGmm}, ...
            {'$s^{{\bf gmm}}$'}, {nGmmPoints}, xylim, cmap);
    end
    
    % Save
    if b_interpEigenvecs
        mVIntToCompare(r,:,:) = VIntToCompare;
        mVNysToCompare(r,:,:) = VNysToCompare;
        mVRefToCompare(r,:,:) = VRefToCompare;
    end
    
end
if b_interpEigenvecs
    vRmseInt = CalcRMSE(mVIntToCompare, mVRefToCompare, 'Analytic');
    vRmseNys = CalcRMSE(mVNysToCompare, mVRefToCompare, 'Nystrom');
    PlotRMSE(sPlotParams, vRmseInt, vRmseNys);
end