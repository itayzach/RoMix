function sWorkspace = Main(presetName, b_saveFigures, b_clearLastRun)
%% Restart
if ~exist('b_clearLastRun', 'var') || b_clearLastRun
    clc; close all; 
end
if ~exist('b_saveFigures', 'var')
    b_saveFigures = false;
end
rng('default'); 
%% GMM / Spectral Clustering
clusterMethod = 'GMM'; % 'GMM' / 'SC'
%% Get perset
sPreset = GetPreset(presetName);
%% Verify preset
VerifyPresetParams(sPreset, clusterMethod);
%% Get plot params
sPlotParams = GetPlotParams(sPreset, b_saveFigures);
%% Run
if sPreset.b_interpEigenvecs
    [tVIntToCompare, tVNysToCompare, tVRepToCompare, tVRefToCompare] = deal(zeros(sPreset.N, sPreset.M, sPreset.R));
    [tVRecPhiToCompare, tVToCompare] = deal(zeros(sPreset.n, sPreset.M, sPreset.R));
end
if sPreset.b_runGraphSignals
    [tSigCnvrtRecPhi, tSigCnvrtRecV, tSigCnvrtRecRep, tSigCnvrt] = deal(zeros(sPreset.n, sPreset.nSignals, sPreset.R));
    [tSigCnvrtInt, tSigCnvrtNys, tSigCnvrtRep, tSigCnvrtRef] = deal(zeros(sPreset.N, sPreset.nSignals, sPreset.R));
end
for r = 1:sPreset.R
    sPlotParams.b_globalPlotEnable = (r == 1) && sPlotParams.b_globalPlotEnable;
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    fprintf('Iteration r = %d of R = %d\n',r,sPreset.R)
    t = tic;
    sDataset = GenerateDataset(sPreset.verticesPDF, sPreset.dim, sPreset.nGenDataCompnts, sPreset.n, sPreset.N, sPreset.dataGenTechnique, sPreset.sDatasetParams);
    xTrain = sDataset.sData.x; yTrain = sDataset.sData.y; xInt = sDataset.sData.xt;
    % ----------------------------------------------------------------------------------------------
    % Estimate distribution parameters
    % ----------------------------------------------------------------------------------------------
    if strcmp(clusterMethod, 'GMM')
        sDistParams = EstimateDistributionParameters(xTrain, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);
        if isfield(sDataset, 'vae')
            sDistParams.vae = sDataset.vae;
        end
    elseif strcmp(clusterMethod, 'SC')
        sDistParams = EstDistParamsSpectClust(xTrain, W, sPreset.gmmNumComponents);
    end
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotClustersAnalysis && sPreset.dim > 1
        PlotGaussianEllipses(sPlotParams, sDistParams);
        PlotCovEigs(sPlotParams, sDistParams);
        if sDistParams.GMModel.NumComponents < 100
            PlotClustersMeans(sPreset, sDistParams);
        end
    end
    % ----------------------------------------------------------------------------------------------
    % Plot dataset
    % ----------------------------------------------------------------------------------------------
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotDataVsGmm 
        nGmmPoints = 50*ismember(sPreset.verticesPDF, {'USPS', 'MNIST'}) + sPreset.n*~ismember(sPreset.verticesPDF, {'USPS', 'MNIST'});
        pltTitle = ['Dataset with n = ', num2str(sPreset.n), ' points'];
        if strcmp(clusterMethod, 'GMM')
            plt2Title = ['Generated ' num2str(nGmmPoints), ' points from ' clusterMethod ' with nEstComp = ' num2str(sPreset.gmmNumComponents)];
        elseif strcmp(clusterMethod, 'SC')
            plt2Title = ['Painted ' num2str(sPreset.n), ' according to ' clusterMethod ' with nEstComp = ' num2str(sPreset.gmmNumComponents)];
        end
        windowStyle = 'normal';
        PlotDataset(sPlotParams, xTrain, yTrain, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle);
    end
    % ----------------------------------------------------------------------------------------------
    % Calculate Phi(xTrain) and lambdaPhi
    % ----------------------------------------------------------------------------------------------
    sKernelParams = CalcKernelParams(sDistParams, sPreset.omegaTilde);
    [sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
         = CalcAnalyticEigenvalues(sPreset.MTilde, sKernelParams);
    [ Phi, lambdaPhi ] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xTrain, sPreset.b_normalizePhi);
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc && sPreset.dim <= 3
        plotInd = [0 4]*(sPreset.dim == 1) + [0 11]*(sPreset.dim > 1);
        figTitle = 'Eigenfunctions of the Gaussian kernel (on $n$ nodes)';
        figName = 'Phi';
        PlotEigenfuncvecScatter([], 'Gaussian', xTrain, [], plotInd(1), plotInd(end), Phi, [], [], [], figTitle, figName, '\phi' );
        figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        PlotSpectrum([], [], lambdaPhi, [], [], '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
    end
    % ----------------------------------------------------------------------------------------------
    % Calculate Phi(xInt)
    % ----------------------------------------------------------------------------------------------
    PhiInt = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xInt, sPreset.b_normalizePhi);
    
    if sPreset.b_compareMethods || sPreset.b_interpEigenvecs || sPlotParams.b_plotMercer
        % ----------------------------------------------------------------------------------------------
        % Build graph - needed for other methods or when interpolating eigenvectors
        % ----------------------------------------------------------------------------------------------
        [WRef, distRef, DRef, LnRef] = CalcAdjacency(xInt, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
        if sPreset.b_takeEigsFromWRef
            W = WRef(1:sPreset.n,1:sPreset.n); dist = distRef(1:sPreset.n,1:sPreset.n); D = DRef(1:sPreset.n,1:sPreset.n); Ln = LnRef(1:sPreset.n,1:sPreset.n);
        else
            [W, dist, D, Ln] = CalcAdjacency(xTrain, sPreset.adjacencyType, sPreset.sDistanceParams, sPreset.omega, sPreset.k, sPreset.nnValue);
        end
        if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotWeights
            PlotWeightsMatrix([], W, dist, D, Ln, xTrain, sPreset.adjacencyType, sPreset.omega, sPreset.k);
        end
        if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotMercer
            CheckMercerTheorem(Phi, lambdaPhi, sPreset.gmmNumComponents, W);
        end
        % ----------------------------------------------------------------------------------------------
        % Perform numeric eigs - needed for other methods or when interpolating eigenvectors
        % ----------------------------------------------------------------------------------------------
        [V, adjLambda, matLambda, VRef, adjLambdaRef, matLambdaRef] = ...
            EigsByTypeWrapper(sPlotParams, sPreset, sDataset, W, D, Ln, WRef, DRef, LnRef);
    else
        Ln = [];
    end
    % ----------------------------------------------------------------------------------------------
    % Eigenvectors interpolation
    % ----------------------------------------------------------------------------------------------
    if sPreset.b_interpEigenvecs
        [tVIntToCompare(:,:,r), tVRecPhiToCompare(:,:,r), tVToCompare(:,:,r), tVRefToCompare(:,:,r), C] = ...
            InterpEigenvectorsRoMix(sPlotParams, sPreset, Phi, lambdaPhi, PhiInt, V, VRef, Ln);

        if sPlotParams.b_globalPlotEnable 
            PlotEigenvectorsWrapper(sPlotParams, sPreset, sDataset, lambdaPhi, Phi, PhiInt, LnRef, tVRefToCompare(:,:,r), tVIntToCompare(:,:,r), C)
        end
        
        if sPreset.b_compareMethods
            [tVNysToCompare(:,:,r)] = InterpEigenvectorsNystrom(sPlotParams, sPreset, sDataset, tVToCompare(:,:,r), adjLambda, tVIntToCompare(:,:,r));
            [tVRepToCompare(:,:,r), mAlpha] = InterpEigenvectorsRepThm(sPlotParams, sPreset, sDataset, W, tVToCompare(:,:,r), tVIntToCompare(:,:,r), Ln);
            if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotC
                CRep = diag(lambdaPhi)*Phi.'*mAlpha;
                PlotCoeffsMatrix(C, '${\bf C}$', CRep, '${\bf C^{\bf rep}} = {\bf \Lambda}{\Phi}^T{\bf \alpha}$', mAlpha, '$\alpha$');
            end
        end

    end
    % ----------------------------------------------------------------------------------------------
    % Graph signals interpolation
    % ----------------------------------------------------------------------------------------------
    if sPreset.b_runGraphSignals
        [tSigCnvrtRecPhi(:,:,r), tSigCnvrt(:,:,r), tSigCnvrtInt(:,:,r), tSigCnvrtRef(:,:,r), mSigCoeffsPhi] = ...
            InterpGraphSignalRoMix(sPreset, sDataset, Phi, lambdaPhi, PhiInt, Ln);

        if sPlotParams.b_globalPlotEnable
            PlotGraphSignalsWrapper(sPlotParams, sPreset, sKernelParams, sDataset, ...
                tSigCnvrt(:,:,r), tSigCnvrtRef(:,:,r), tSigCnvrtRecPhi(:,:,r), tSigCnvrtInt(:,:,r), mSigCoeffsPhi)
        end

        if sPreset.b_compareMethods
            [tSigCnvrtRecRep(:,:,r), tSigCnvrtRep(:,:,r)] = InterpGraphSignalRepThm(sPreset, sDataset, W, Ln);
            [tSigCnvrtRecV(:,:,r), tSigCnvrtNys(:,:,r)] = InterpGraphSignalNystrom(sPreset, sDataset, V, adjLambda);
        end
    end
    fprintf('Iteration r = %d of R = %d finished (took %.2f sec)\n',r,sPreset.R,toc(t))
end
% ------------------------------------------------------------------------------------------
% Accuracy
% ------------------------------------------------------------------------------------------
if sPreset.b_interpEigenvecs
    InterpEigenvecsMetrics(sPlotParams, sPreset, tVIntToCompare, tVRefToCompare, tVRecPhiToCompare, tVToCompare, tVNysToCompare, tVRepToCompare);
end
if sPreset.b_runGraphSignals
    InterpGraphSignalsMetrics(sPlotParams, sPreset, tSigCnvrtRecPhi, tSigCnvrt, tSigCnvrtInt, tSigCnvrtRef, tSigCnvrtRecRep, tSigCnvrtRep, tSigCnvrtRecV, tSigCnvrtNys);
end
%% Save workspace
sWorkspace = SaveWorkspaceToStruct();
ReadWorkspaceStructToBase(sWorkspace);
end