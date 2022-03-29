clear; clc;
close all;
M = 20;

%% 1D Grid analytic
Nx1d = 2500;
xMax1d = 5;
len = xMax1d;
grid1dData = GenerateGridData(numel(Nx1d), Nx1d, 0, xMax1d);
[grid1dV, ~, grid1dMatLambda] = CalcAnalyticLapEigsGrid(grid1dData, M, len);
PlotEigenfuncvecScatter([], 'Grid', grid1dData, [], 0, 4, ...
   grid1dV, grid1dMatLambda, '\lambda^{v}', [], 'Grid - Eigenvectors', 'figName', 'v');
PlotInnerProductMatrix([], grid1dV, [], 'Grid 1d - ${\bf V}^T {\bf V}$', 'V');

%% Swiss Roll
N = 1000;
b_plotSwissRoll = false;
height = 20;
maxTheta = 4*pi;

[srData, S] = GenerateSwissRoll(N);
len = [SwissRollArclength(maxTheta), height];
[srV, srAdjLambda, srMatLambda] = CalcAnalyticLapEigsGrid(S, M, len);
PlotEigenfuncvecScatter([], 'SwissRoll', S, [], 0, 8, ...
   srV, srMatLambda, '\lambda^{v}', [], 'Swiss Roll - Parametrization', 'figName', 'v');
PlotEigenfuncvecScatter([], 'SwissRoll', srData, [], 0, 8, ...
   srV, srMatLambda, '\lambda^{v}', [], 'Swiss Roll - Eigenvectors', 'figName', 'v');
PlotInnerProductMatrix([], srV, [], 'Swiss Roll - ${\bf V}^T {\bf V}$', 'V');

%% 2D Grid analytic
Nx2d = [10 10];
xMax2d = [1 1];

grid2dData = GenerateGridData(numel(Nx2d), Nx2d, [0 0], xMax2d);
len = xMax2d;
[grid2dV, ~, grid2dMatLambda] = CalcAnalyticLapEigsGrid(grid2dData, M, len);
PlotEigenfuncvecScatter([], 'Grid', grid2dData, [], 0, 8, ...
   grid2dV, grid2dMatLambda, '\lambda^{v}', [], 'Grid - Eigenvectors', 'figName', 'v');
PlotInnerProductMatrix([], grid2dV, [], 'Grid 2d - ${\bf V}^T {\bf V}$', 'V');



%% 3D Grid analytic
Nx3d = [20 20 20];
xMax3d = [1 1 1];
len = xMax3d;
grid3dData = GenerateGridData(numel(Nx3d), Nx3d, [0 0 0], xMax3d);
[grid3dV, ~, grid3dMatLambda] = CalcAnalyticLapEigsGrid(grid3dData, M, len);
PlotInnerProductMatrix([], grid3dV, [], 'Grid 3d - ${\bf V}^T {\bf V}$', 'V');
