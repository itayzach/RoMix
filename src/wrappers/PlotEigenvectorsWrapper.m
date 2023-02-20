function PlotEigenvectorsWrapper(sPlotParams, sPreset, sDataset, mSigCnvrtRecRef, mSigCnvrtIntRef, mSigCnvrtRec, mSigCnvrtInt)
xTrain = sDataset.sData.x;
xInt = sDataset.sData.xt;
% ----------------------------------------------------------------------------------------------
% Plots
% ----------------------------------------------------------------------------------------------
if sPlotParams.b_globalPlotEnable && sPreset.dim <= 3
    if sPreset.dim == 1
        plotInd = [0:3, sPreset.M-1];
        cMarkers = cell(1,2*numel(plotInd));
        [cMarkers{1:numel(plotInd)}] = deal('o');
        [cMarkers{numel(plotInd)+1:numel(plotInd)*2}] = deal('.');
    else
        plotInd = round(linspace(1,sPreset.M-1,4));
        [cMarkers{1:numel(plotInd)*2}] = deal('.');
    end

    % Interpolation
    [cData{1:numel(plotInd)*2}] = deal(xTrain);
    cSigStr = [RepLegend('\\psi', plotInd), RepLegend('\\psi^{{\\bf RoMix}}', plotInd)];
    [cNumCircles{1:numel(plotInd)*2}] = deal((1:sPreset.nLabeled).');
    pltTitle = []; %['Eigenvectors interpolation ($\ell = ', num2str(sPreset.nLabeled), '$, $n = ', num2str(sPreset.n), '$)'];
    PlotGraphSignals(sPlotParams, pltTitle, ...
        [sPreset.matrixForEigs, '_InterpEigs_',num2str(plotInd(1)), '_to_', num2str(plotInd(end)) ], cData, ...
        [mat2cell(mSigCnvrtRecRef(:,plotInd+1),sPreset.n,ones(1,numel(plotInd))), mat2cell(mSigCnvrtRec(:,plotInd+1),sPreset.n,ones(1,numel(plotInd)))], ...
        cSigStr, cNumCircles, cMarkers, [], [], [min(min(mSigCnvrtRecRef(:,plotInd+1))), max(max(mSigCnvrtRecRef(:,plotInd+1)))]);

    % Extrapolation
    [cData{1:numel(plotInd)*2}] = deal(xInt);
    cSigStr = [RepLegend('\\tilde{\\psi}', plotInd), RepLegend('\\tilde{\\psi}^{{\\bf RoMix}}', plotInd)];
    [cNumCircles{1:numel(plotInd)*2}] = deal((1:sPreset.N).');
    pltTitle = []; %['Eigenvectors extrapolation ($N = ', num2str(sPreset.N), '$)'];
    PlotGraphSignals(sPlotParams, pltTitle, ...
        [sPreset.matrixForEigs, '_ExtrapEigs_',num2str(plotInd(1)), '_to_', num2str(plotInd(end)) ], cData, ...
        [mat2cell(mSigCnvrtIntRef(:,plotInd+1),sPreset.N,ones(1,numel(plotInd))), mat2cell(mSigCnvrtInt(:,plotInd+1),sPreset.N,ones(1,numel(plotInd)))], ...
        cSigStr, cNumCircles, cMarkers, [], [], [min(min(mSigCnvrtIntRef(:,plotInd+1))), max(max(mSigCnvrtIntRef(:,plotInd+1)))]);
end
end