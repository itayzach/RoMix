function sSimParams = GetSimParams(M)
%% Num of eigenfunctions
if ~exist('M', 'var')
    sSimParams.CalcEigenFuncsM = 12;
    sSimParams.PlotEigenFuncsM = 12;
    sSimParams.PlotSpectM = 12;
else
    sSimParams.CalcEigenFuncsM = M;
    sSimParams.PlotEigenFuncsM = M;
    sSimParams.PlotSpectM = M;
end
%% Figures
sSimParams.outputFolder = 'figs';
sSimParams.b_plotEigenfunctions = true;
end