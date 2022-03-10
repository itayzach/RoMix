function sWorkspace = Main(presetName, b_clearLastRun)
%% Restart
if ~exist('b_clearLastRun', 'var') || b_clearLastRun
    clc; close all; set(0,'DefaultFigureWindowStyle','normal')
end
rng('default'); 
%% Get perset from input
if exist('presetName', 'var')
    if isstruct(presetName)
        sPreset = presetName;
    elseif ischar(presetName)
        sPreset = eval(presetName);
    end
end
PrintPreset(sPreset);
%% GMM / Spectral Clustering
clusterMethod = 'GMM'; % 'GMM' / 'SC'
%% Verify
assert(~sPreset.b_debugUseAnalytic || ...
    (sPreset.b_debugUseAnalytic && ismember(sPreset.verticesPDF, {'Gaussian', 'Grid', 'Uniform', 'SwissRoll'})))
assert(~strcmp(sPreset.adjacencyType,'NearestNeighbor') || ...
    strcmp(sPreset.adjacencyType,'NearestNeighbor') && strcmp(sPreset.verticesPDF,'Grid'))
assert((strcmp(clusterMethod, 'GMM') && sPreset.n >= sPreset.dim) || (strcmp(clusterMethod, 'SC')))
%% Plot params
sPlotParams = GetPlotParams(sPreset);
%% Run
[tVIntToCompare, tVNysToCompare, tVRepToCompare, tVRefToCompare] = deal(zeros(sPreset.N, sPreset.M, sPreset.R));
[tSigCnvrtRecPhi, tSigCnvrtRecV, tSigCnvrtRecRep, tSigCnvrt] = deal(zeros(sPreset.n, sPreset.nSignals, sPreset.R));
[tSigCnvrtInt, tSigCnvrtNys, tSigCnvrtRep, tSigCnvrtRef] = deal(zeros(sPreset.N, sPreset.nSignals, sPreset.R));
for r = 1:sPreset.R
    sPlotParams.b_globalPlotEnable = (r == 1) && sPlotParams.b_globalPlotEnable;
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    fprintf('Iteration r = %d of R = %d\n',r,sPreset.R)
    sDataset = GenerateDataset(sPreset.verticesPDF, sPreset.dim, sPreset.nGenDataCompnts, sPreset.n, sPreset.N, sPreset.dataGenTechnique, sPreset.sDatasetParams);
    xTrain = sDataset.sData.x; yTrain = sDataset.sData.y; xInt = sDataset.sData.xt;
    % ----------------------------------------------------------------------------------------------
    % Build graph (not needed for EigsRLS)
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
    % ----------------------------------------------------------------------------------------------
    % Estimate distribution parameters
    % ----------------------------------------------------------------------------------------------
    if strcmp(clusterMethod, 'GMM')
        sDistParams = EstimateDistributionParameters(xTrain, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);
    elseif strcmp(clusterMethod, 'SC')
        sDistParams = EstDistParamsSpectClust(xTrain, W, sPreset.gmmNumComponents);
    end
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotClustersAnalysis
        PlotGaussianEllipses(sPlotParams, sDistParams);
        PlotCovEigs(sPlotParams, sDistParams);
        PlotClustersMeans(sPreset, sDistParams);
    end
    % ----------------------------------------------------------------------------------------------
    % Plot dataset
    % ----------------------------------------------------------------------------------------------
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotDataVsGmm && sPreset.dim <= 3
        nGmmPoints = sPreset.n;
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
    % Perform numeric eigs
    % ----------------------------------------------------------------------------------------------
    [V, adjLambda, matLambda, VRef, adjLambdaRef, matLambdaRef] = ...
        EigsByTypeWrapper(sPlotParams, sPreset, sDataset, W, D, Ln, WRef, DRef, LnRef);
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
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotMercer
        CheckMercerTheorem(Phi, lambdaPhi, sPreset.gmmNumComponents, W);
    end
    % ----------------------------------------------------------------------------------------------
    % Calculate Phi(xInt)
    % ----------------------------------------------------------------------------------------------
    PhiInt = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xInt, sPreset.b_normalizePhi);
    % ----------------------------------------------------------------------------------------------
    % Calculate Nystrom eigenvectors
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = CalcDistance(xTrain, xInt, sPreset.sDistanceParams);
    WTrainInt = exp(-distLUBlockRUBlock.^2/(2*sPreset.omega^2));
    interpRatio = sPreset.N/sPreset.n;
    lambdaNys = adjLambda*sqrt(interpRatio);
    VNys = WTrainInt.'*V*diag(1./lambdaNys);
    % ----------------------------------------------------------------------------------------------
    % Eigenvectors interpolation
    % ----------------------------------------------------------------------------------------------
    if sPreset.b_interpEigenvecs
        [tVIntToCompare(:,:,r), tVNysToCompare(:,:,r), tVRepToCompare(:,:,r), tVRefToCompare(:,:,r)] = ...
            InterpEigenvectors(sPlotParams, sPreset, sDataset, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, VRef, W, WRef, D, DRef, Ln, LnRef);
    end
    % ----------------------------------------------------------------------------------------------
    % Graph signals interpolation
    % ----------------------------------------------------------------------------------------------
    if sPreset.b_runGraphSignals
        [tSigCnvrtRecPhi(:,:,r), tSigCnvrtRecV(:,:,r), tSigCnvrtRecRep(:,:,r), tSigCnvrt(:,:,r), ...
         tSigCnvrtInt(:,:,r), tSigCnvrtNys(:,:,r), tSigCnvrtRep(:,:,r), tSigCnvrtRef(:,:,r)] = ...
            InterpGraphSignal(sPlotParams, sPreset, sDataset, sKernelParams, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, W, WRef, D, DRef, Ln, LnRef);
    end
end
% ------------------------------------------------------------------------------------------
% Accuracy
% ------------------------------------------------------------------------------------------
if sPreset.b_interpEigenvecs
    [vRmseInt, vMseInt, vAccInt, mErrInt, vCohInt, vAccStdInt] = CalcErrAndAcc(tVIntToCompare, tVRefToCompare, 'Analytic');
    if sPreset.b_compareMethods
        [vRmseNys, vMseNys, vAccNys, mErrNys, vCohNys, vAccStdNys] = CalcErrAndAcc(tVNysToCompare, tVRefToCompare, 'Nystrom');
        [vRmseRep, vMseRep, vAccRep, mErrRep, vCohRep, vAccStdRep] = CalcErrAndAcc(tVRepToCompare, tVRefToCompare, 'Representer');
        PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], [vAccStdInt, vAccStdNys, vAccStdRep], ...
            {'Acc$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'Acc$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
             'Acc$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'}, ['Acc_eigs_0_to_' num2str(sPreset.M-1)]);
    else
        PlotAccuracy(sPlotParams, vAccInt, vAccStdInt, {'Acc$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$'}, ['Acc_eigs_0_to_' num2str(sPreset.M-1)]);
    end
end
if sPreset.b_runGraphSignals
    [vRmseRecPhi, vMseRecPhi, vAccRecPhi, mErrRecPhi, vCohRecPhi, vAccStdPhi] = CalcErrAndAcc(tSigCnvrtRecPhi, tSigCnvrt, 'EigsRLS (train)');
    [vRmseInt, vMseInt, vAccInt, mErrInt, vCohInt, vAccStdInt]                = CalcErrAndAcc(tSigCnvrtInt, tSigCnvrtRef, 'EigsRLS (test)');
    if sPreset.b_compareMethods
        [vRmseRecRep, vMseRecRep, vAccRecRep, mErrRecRep, vCohRecRep, vAccStdRecRep] = CalcErrAndAcc(tSigCnvrtRecRep, tSigCnvrt, 'Representer (train)');
        [vRmseRep, vMseRep, vAccRep, mErrRep, vCohRep, vAccStdRep]                   = CalcErrAndAcc(tSigCnvrtRep, tSigCnvrtRef, 'Representer (test)');
        [vRmseRecV, vMseRecV, vAccRecV, mErrRecV, vCohRecV, vAccStdRecV]             = CalcErrAndAcc(tSigCnvrtRecV, tSigCnvrt, 'Nystrom (train)');
        [vRmseNys, vMseNys, vAccNys, mErrNys, vCohNys, vAccStdNys]                   = CalcErrAndAcc(tSigCnvrtNys, tSigCnvrtRef, 'Nystrom (test)');
        if isfield(sPreset.sDatasetParams, 'monthNames')
            PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], [vAccStdInt, vAccStdNys, vAccStdRep], ...
                {'Acc$(s^{{\bf int}}_m, s^{{\bf ref}}_m)$', 'Acc$(s^{{\bf nys}}_m, s^{{\bf ref}}_m)$', ...
                'Acc$(s^{{\bf rep}}_m, s^{{\bf ref}}_m)$'}, 'InterpAcc', [], sPreset.sDatasetParams.monthNames, 'Interpolation accuracy');
            PlotAccuracy(sPlotParams, [vAccRecPhi, vAccRecV, vAccRecRep], [vAccStdPhi, vAccStdRecV, vAccStdRecRep], ...
                {'Acc$(s^{{\bf int}}_m, s^{{\bf ref}}_m)$', 'Acc$(s^{{\bf nys}}_m, s^{{\bf ref}}_m)$', ...
                'Acc$(s^{{\bf rep}}_m, s^{{\bf ref}}_m)$'}, 'ProjAcc', [], sPreset.sDatasetParams.monthNames, 'Projection accuracy');
        end
    end
end
sWorkspace = SaveWorkspaceToStruct();
ReadWorkspaceStructToBase(sWorkspace);
end