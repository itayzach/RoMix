function [sParams] = GetOldParameters()

% %% x-axis
% sParams.n_x_axis = 100;
% sParams.dx = (sParams.xMax - sParams.xMin)/sParams.n_x_axis;
% x = zeros(sParams.n_x_axis, sParams.dim);
% for d = 1:sParams.dim
%     x(:,d) = (sParams.xMin(d) : sParams.dx(d) : sParams.xMax(d)-sParams.dx(d)).';
% end
% sParams.x = x;


%% num of eigenfunctions
sParams.PlotEigenFuncsM = 20;
sParams.PlotSpectM = 50;
sParams.RkhsM = 20;
sParams.OrthM = 30;
sParams.MercerM = 0;
sParams.M = 20;
sParams.FirstM = 0;

%% simulation
sParams.sSim.outputFolder = 'figs';

sParams.sSim.b_plotEigenFigs          = true;
sParams.sSim.b_verifyKernelEigenfuncs = false;
sParams.sSim.b_verifyEigOrth          = false;
sParams.sSim.b_verifyMercersTheorem   = false;
sParams.sSim.b_extrapolateEnable      = false;

sParams.sSim.b_randomStepSize = true;
sParams.sSim.b_plot_contourf = false;

end