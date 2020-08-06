function sSimParams = GetSimParams()

%% Num of eigenfunctions
sSimParams.CalcEigenFuncsM = 50;
sSimParams.PlotEigenFuncsM = 20;
sSimParams.PlotSpectM = 50;

%% Figures
sSimParams.outputFolder = 'figs';
sSimParams.b_plotEigenfunctions = true;
end