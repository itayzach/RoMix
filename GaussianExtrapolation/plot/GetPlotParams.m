function sPlotParams = GetPlotParams()
%% Save figs
sPlotParams.outputFolder = 'figs';
sPlotParams.b_globalPlotEnable          = true;
%% Eigenvectors
sPlotParams.b_plotOrigEvecs             = false;
sPlotParams.b_plotGaussianKernelEigfunc = false;
%% Debug
sPlotParams.b_plotWeights               = false;
sPlotParams.b_plotC                     = false;
sPlotParams.b_plotInnerProductMatrices  = false;
sPlotParams.b_plotMercer                = false;
%% Data
sPlotParams.b_plotDataVsGmm             = true;
%% Graph signals
sPlotParams.b_plotGmmSignal             = false;
sPlotParams.b_plotExtraGraphSigAnalysis = false;

end