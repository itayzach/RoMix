%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
outputFolder = 'figs'; figSaveType = 'png';

%% Random data
verticesPDF = 'Uniform_1D'; % 'Gaussian_1D' / 'Uniform_1D'
n = 1000;
if strcmp(verticesPDF, 'Uniform_1D')
    xMax = 10;
    xMin = -10;
    xTrain = (xMax - xMin)*rand(n,1) + xMin; %linspace(minVal,maxVal,N)';
elseif strcmp(verticesPDF, 'Gaussian_1D')
    sigma = 10;
    mu = 0;
    xMax = 3*sigma;
    xMin = -3*sigma;
    xTrain = sigma*randn(n,1) + mu;
else
    error('invalid verticesPDF');
end

fig = figure('Name', 'Histogram of X'); 
histogram(xTrain,100);
title('Histogram of ${\bf X}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig1_histogram_X'), figSaveType);
%% Generate graph
omega = 3;
dist = pdist2(xTrain, xTrain);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
M = 20;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:5;

fig = figure('Name', 'Eigenvectors of W');
plot(xTrain, V(:,vInd),'.');
title('Eigenvectors of ${\bf W}$ ', 'interpreter', 'latex', 'FontSize', 16);
legend(strcat('$v_{',string(vInd),'}$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig2_evecs_W'), figSaveType);

%% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
[estCdf_xTrainGrid, xTrainGrid] = ecdf(xTrain);
xTrainGrid = xTrainGrid(1:end-1);
estCdf_xTrainGrid = estCdf_xTrainGrid(1:end-1);

pCdfDegree = 10;
invpCdfDegree = 10;
pCdf = polyfit(xTrainGrid, estCdf_xTrainGrid, pCdfDegree); % to have analytic expression for the cdf
invpCdf = polyfit(estCdf_xTrainGrid, xTrainGrid, invpCdfDegree); % to have analytic expression for the cdf

%% Transform to Gtilde
muTilde = 0;
sigmaTilde = 1;
xTildeTrain = T(pCdf, true, muTilde, sigmaTilde, xTrain);

%% Demonstrate T (1/2)
PlotPolyCdfDemonstration1(xMin, xMax, pCdf, xTrainGrid, estCdf_xTrainGrid, muTilde, sigmaTilde)

%% Demonstrate T (2/2)
PlotPolyCdfDemonstration2(xMin, xMax, pCdf, invpCdf, muTilde, sigmaTilde)

%% Build G tilde
omegaTilde = 0.3;
distTilde = pdist2(xTildeTrain, xTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));

%% Calculate analytic eigenfunctions of W_tilde, and numeric eigenvectors for comparison
MTilde = 20;
[PhiTilde, ~] = SimpleCalcAnalyticEigenfunctions(xTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

fig = figure('Name', 'evecs vs. efuncs of WTilde');
plot(xTildeTrain, VTilde(:,vInd),'o');
hold on
plot(xTildeTrain, PhiTilde(:,vInd),'.');
legend([strcat('$\tilde{v}_{',string(vInd),'}$') strcat('$\tilde{\phi}_',string(vInd),'$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('Eigenfunctions vs. eigenvectors of $\tilde{{\bf W}}$ (from $x_{{\bf train}})$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig5_evecs_vs_efuncs'), figSaveType);

%% Functional maps
C = pinv(PhiTilde)*V;
pinvC = pinv(C);

fig = figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C} = \tilde{{\bf \Phi}}^\dagger {\bf V}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig6_C'), figSaveType);

fig = figure('Name', 'pinv(C)');
imagesc(pinvC);
colorbar();
title('${\bf C}^\dagger$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig7_pinvC'), figSaveType);

%% V in terms of Phi_tilde
xEstTrain = invT(invpCdf, muTilde, sigmaTilde, xTildeTrain);
VRec = PhiTilde*C;

fig = figure('Name', 'V in terms of PhiTilde');
plot(xTrain, V(:,vInd),'o');
hold on
plot(xEstTrain, VRec(:,vInd),'.');
legend([strcat('$v_{',string(vInd),'}$') strcat('$v_{{\bf rec},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['${\bf V}$ in terms of ${\bf \tilde{\Phi}}$' newline '${\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig8_VRec'), figSaveType);

%% Interpolate
N = 5000;
xInt = linspace(xMin+0.5, xMax-0.5, N)';
xTildeInt = T(pCdf, true, muTilde, sigmaTilde, xInt);
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(xTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);

fig = figure('Name', 'Eigenfunctions of WTilde on the entire axis');
plot(xTildeInt, PhiTildeInt(:,vInd),'LineWidth',2);
legend(strcat('$\tilde{\phi}_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Eigenfunctions of $\tilde{{\bf W}}$ on the entire axis' newline ...
    '$\tilde{\omega}$ = ' num2str(omegaTilde, '%.2f') ...
    '; $\tilde{\sigma}$ = ' num2str(sigmaTilde, '%.2f') ...
    '; $\tilde{\mu}$ = ' num2str(muTilde, '%.2f') ], 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig9_efuncs_WTilde_int'), figSaveType);

xIntInvT = invT(invpCdf, muTilde, sigmaTilde, xTildeInt);
VInt = PhiTildeInt*C;

interpRatio = N/n;
VIntRenormed = sqrt(interpRatio)*VInt;

fig = figure('Name', 'Interpolated evecs of W');
plot(xTrain, V(:,vInd),'o');
hold on
plot(xIntInvT, VIntRenormed(:,vInd),'.');
legend([strcat('$v_{',string(vInd),'}$') strcat('$v_{{\bf int},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Interpolated eigenvectors of ${\bf W}$' newline '${\bf V}_{{\bf int}} = \sqrt{\frac{N}{n}}{\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig10_VInt'), figSaveType);
