function Main(presetName)
%% Restart
% clear;
clc; rng('default'); 
% close all; 
set(0,'DefaultFigureWindowStyle','normal')
%% Load preset
% sPreset = Get1DGridPreset();
% sPreset = Get2DGridPreset();
% sPreset = Get1DUniformPreset();
% sPreset = Get1DGaussianPreset();
% sPreset = GetSwissRollPreset(); % For graph sig interp, run with n=1000
% sPreset = GetBrazilWeatherPreset();
% sPreset = GetTwoMoonsPreset(); %Classifier is okay. Should fix order of VRef relative to V...
% sPreset = GetUspsPreset();
% sPreset = GetMnistPreset();
%% Get perset from input
if exist('presetName', 'var')
    if isstruct(presetName)
        sPreset = presetName;
    elseif ischar(presetName)
        sPreset = eval(presetName);
    end
end
PrintPresetName(sPreset);
%%
clusterMethod = 'GMM'; % 'GMM' / 'SC'
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
nSig               = sPreset.nSignals;
R                  = sPreset.R;
b_debugUseAnalytic = sPreset.b_debugUseAnalytic;
b_forceCtoIdentity = sPreset.b_forceCtoIdentity;
b_normalizePhi     = sPreset.b_normalizePhi;
b_takeEigsFromWRef = sPreset.b_takeEigsFromWRef;
b_flipSign         = sPreset.b_flipSign;
b_pairwiseFlipSign = sPreset.b_pairwiseFlipSign;
b_interpEigenvecs  = sPreset.b_interpEigenvecs;
b_runGraphSignals  = sPreset.b_runGraphSignals;
b_compareMethods   = sPreset.b_compareMethods;
interpRatio        = N/n;
%% Verify
assert(~b_debugUseAnalytic || ...
    (b_debugUseAnalytic && ismember(verticesPDF, {'Gaussian', 'Grid', 'Uniform', 'SwissRoll'})))
assert(~strcmp(adjacencyType,'NearestNeighbor') || ...
    strcmp(adjacencyType,'NearestNeighbor') && strcmp(verticesPDF,'Grid'))
assert((strcmp(clusterMethod, 'GMM') && n >= dim) || (strcmp(clusterMethod, 'SC')))
%% Plot params
sPlotParams = GetPlotParams();
sPlotParams.actualDataDist = verticesPDF;
sPlotParams.matrixForEigs = matrixForEigs;
sPlotParams.dim = dim;
%% Run
[tVIntToCompare, tVNysToCompare, tVRepToCompare, tVRefToCompare] = deal(zeros(N, M, R));
[tSigCnvrtRecPhi, tSigCnvrtRecV, tSigCnvrtRecRep, tSigCnvrt] = deal(zeros(n, nSig, R));
[tSigCnvrtInt, tSigCnvrtNys, tSigCnvrtRep, tSigCnvrtRef] = deal(zeros(N, nSig, R));
for r = 1:R
    sPlotParams.b_globalPlotEnable = (r == 1) && sPlotParams.b_globalPlotEnable;
    % ----------------------------------------------------------------------------------------------
    % Generate dataset
    % ----------------------------------------------------------------------------------------------
    fprintf('Iteration r = %d of R = %d\n',r,R)
    sDataset = GenerateDataset(verticesPDF, dim, nGenDataCompnts, n, N, dataGenTechnique, sDatasetParams);
    xTrain = sDataset.sData.x; yTrain = sDataset.sData.y; xInt = sDataset.sData.xt;
    % ----------------------------------------------------------------------------------------------
    % Build graph (not needed for EigsRLS)
    % ----------------------------------------------------------------------------------------------
    [WRef, distRef, DRef, LnRef] = SimpleCalcAdjacency(xInt, adjacencyType, sDistanceParams, omega, k, nnValue);
    if b_takeEigsFromWRef
        W = WRef(1:n,1:n); dist = distRef(1:n,1:n); D = DRef(1:n,1:n); Ln = LnRef(1:n,1:n);
    else
        [W, dist, D, Ln] = SimpleCalcAdjacency(xTrain, adjacencyType, sDistanceParams, omega, k, nnValue);
    end
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotWeights
        PlotWeightsMatrix([], W, dist, D, Ln, xTrain, adjacencyType, omega, k);
    end
    % ----------------------------------------------------------------------------------------------
    % Estimate distribution parameters
    % ----------------------------------------------------------------------------------------------
    if strcmp(clusterMethod, 'GMM')
        sDistParams = EstimateDistributionParameters(xTrain, gmmNumComponents, gmmRegVal, gmmMaxIter);
    elseif strcmp(clusterMethod, 'SC')
        sDistParams = EstDistParamsSpectClust(xTrain, W, gmmNumComponents);
    end
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotDataVsGmm && dim <= 3
        nGmmPoints = n;
        pltTitle = ['Dataset with n = ', num2str(n), ' points'];
        if strcmp(clusterMethod, 'GMM')
            plt2Title = ['Generated ' num2str(nGmmPoints), ' points from ' clusterMethod ' with nEstComp = ' num2str(gmmNumComponents)];
        elseif strcmp(clusterMethod, 'SC')
            plt2Title = ['Painted ' num2str(n), ' according to ' clusterMethod ' with nEstComp = ' num2str(gmmNumComponents)];
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
    sKernelParams = CalcKernelParams(sDistParams, omegaTilde);
    [sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
         = CalcAnalyticEigenvalues(MTilde, sKernelParams);
    [ Phi, lambdaPhi ] = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xTrain, b_normalizePhi);
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotGaussianKernelEigfunc && dim <= 3
        figTitle = 'Eigenfunctions of the Gaussian kernel (on $n$ nodes)';
        figName = 'Phi';
        PlotEigenfuncvecScatter([], 'Gaussian', xTrain, [], 0, 4, Phi, [], [], [], figTitle, figName, '\phi' );
        figTitle = 'Analytic eigenvalues of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$';
        PlotSpectrum([], [], lambdaPhi, [], [], '\tilde{\lambda}^{\phi}_m', [], [], figTitle);
    end
    if sPlotParams.b_globalPlotEnable && sPlotParams.b_plotMercer
        CheckMercerTheorem(Phi, lambdaPhi, gmmNumComponents, W);
    end
    % ----------------------------------------------------------------------------------------------
    % Calculate Phi(xInt)
    % ----------------------------------------------------------------------------------------------
    PhiInt = CalcAnalyticEigenfunctions(MTilde, sKernelParams, xInt, b_normalizePhi);
    % ----------------------------------------------------------------------------------------------
    % Calculate Nystrom eigenvectors
    % ----------------------------------------------------------------------------------------------
    distLUBlockRUBlock = CalcDistance(xTrain, xInt, sDistanceParams);
    WTrainInt = exp(-distLUBlockRUBlock.^2/(2*omegaNys^2));
    lambdaNys = adjLambda*sqrt(interpRatio);
    VNys = WTrainInt.'*V*diag(1./lambdaNys);
    % ----------------------------------------------------------------------------------------------
    % Eigenvectors interpolation
    % ----------------------------------------------------------------------------------------------
    if b_interpEigenvecs
        [tVIntToCompare(:,:,r), tVNysToCompare(:,:,r), tVRepToCompare(:,:,r), tVRefToCompare(:,:,r)] = ...
            InterpEigenvectors(sPlotParams, sPreset, sDataset, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, VRef, W, WRef, D, DRef, Ln, LnRef);
    end
    % ----------------------------------------------------------------------------------------------
    % Graph signals interpolation
    % ----------------------------------------------------------------------------------------------
    if b_runGraphSignals
        [tSigCnvrtRecPhi(:,:,r), tSigCnvrtRecV(:,:,r), tSigCnvrtRecRep(:,:,r), tSigCnvrt(:,:,r), ...
         tSigCnvrtInt(:,:,r), tSigCnvrtNys(:,:,r), tSigCnvrtRep(:,:,r), tSigCnvrtRef(:,:,r)] = ...
            InterpGraphSignal(sPlotParams, sPreset, sDataset, sKernelParams, Phi, lambdaPhi, PhiInt, VNys, WTrainInt, V, W, WRef, D, DRef, Ln, LnRef);
    end
end
% ------------------------------------------------------------------------------------------
% Accuracy
% ------------------------------------------------------------------------------------------
if b_interpEigenvecs
    [vRmseInt, vMseInt, vAccInt, mErrInt, vCohInt, vAccStdInt] = CalcErrAndAcc(tVIntToCompare, tVRefToCompare, 'Analytic');
    if b_compareMethods
        [vRmseNys, vMseNys, vAccNys, mErrNys, vCohNys, vAccStdNys] = CalcErrAndAcc(tVNysToCompare, tVRefToCompare, 'Nystrom');
        [vRmseRep, vMseRep, vAccRep, mErrRep, vCohRep, vAccStdRep] = CalcErrAndAcc(tVRepToCompare, tVRefToCompare, 'Representer');
        PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], [vAccStdInt, vAccStdNys, vAccStdRep], ...
            {'Acc$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$', 'Acc$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$', ...
             'Acc$(v^{{\bf rep}}_m, v^{{\bf ref}}_m)$'}, ['Acc_eigs_0_to_' num2str(M-1)]);
    else
        PlotAccuracy(sPlotParams, vAccInt, vAccStdInt, {'Acc$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$'}, ['Acc_eigs_0_to_' num2str(M-1)]);
    end
end
if b_runGraphSignals
    [vRmseRecPhi, vMseRecPhi, vAccRecPhi, mErrRecPhi, vCohRecPhi, vAccStdPhi] = CalcErrAndAcc(tSigCnvrtRecPhi, tSigCnvrt, 'EigsRLS (train)');
    [vRmseInt, vMseInt, vAccInt, mErrInt, vCohInt, vAccStdInt]                = CalcErrAndAcc(tSigCnvrtInt, tSigCnvrtRef, 'EigsRLS (test)');
    if b_compareMethods
        [vRmseRecRep, vMseRecRep, vAccRecRep, mErrRecRep, vCohRecRep, vAccStdRecRep] = CalcErrAndAcc(tSigCnvrtRecRep, tSigCnvrt, 'Representer (train)');
        [vRmseRep, vMseRep, vAccRep, mErrRep, vCohRep, vAccStdRep]                   = CalcErrAndAcc(tSigCnvrtRep, tSigCnvrtRef, 'Representer (test)');
        [vRmseRecV, vMseRecV, vAccRecV, mErrRecV, vCohRecV, vAccStdRecV]             = CalcErrAndAcc(tSigCnvrtRecV, tSigCnvrt, 'Nystrom (train)');
        [vRmseNys, vMseNys, vAccNys, mErrNys, vCohNys, vAccStdNys]                   = CalcErrAndAcc(tSigCnvrtNys, tSigCnvrtRef, 'Nystrom (test)');
        if isfield(sDatasetParams, 'monthNames')
            PlotAccuracy(sPlotParams, [vAccInt, vAccNys, vAccRep], [vAccStdInt, vAccStdNys, vAccStdRep], ...
                {'Acc$(s^{{\bf int}}_m, s^{{\bf ref}}_m)$', 'Acc$(s^{{\bf nys}}_m, s^{{\bf ref}}_m)$', ...
                'Acc$(s^{{\bf rep}}_m, s^{{\bf ref}}_m)$'}, 'InterpAcc', [], sDatasetParams.monthNames, 'Interpolation accuracy');
            PlotAccuracy(sPlotParams, [vAccRecPhi, vAccRecV, vAccRecRep], [vAccStdPhi, vAccStdRecV, vAccStdRecRep], ...
                {'Acc$(s^{{\bf int}}_m, s^{{\bf ref}}_m)$', 'Acc$(s^{{\bf nys}}_m, s^{{\bf ref}}_m)$', ...
                'Acc$(s^{{\bf rep}}_m, s^{{\bf ref}}_m)$'}, 'ProjAcc', [], sDatasetParams.monthNames, 'Projection accuracy');
        end
    end
end

end