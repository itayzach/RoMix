%% Restart
clc; clear; close all; 
rng('default')

%% Parameters
graphName =  'Uniform_1D'; % Uniform_1D
N = 1000; % number of nodes in G

%==========================================================================
% kernel parameters
%==========================================================================
omega       = 0.3;
estDataDist = 'Gaussian';
nComponents = 1;

%==========================================================================
% basis parameters
%==========================================================================
M_G = 5;
M_Gtilde = 30;

vInd = 4;     % one-based
v_m = vInd-1; % zero-based

C_TransformType   = 'PhiTilde^(-1)*R*Phi'; % 'PhiTilde^(-1)*R*Phi' / 'PhiTilde^(-1)*Phi'
R_TransformType   = 'permutation'; % 'permutation' / 'randomMatrix' / 'eye'
phiTildeType      = 'Analytic'; % 'Analytic' / 'Numeric' / 'Compare'
b_runInterpolation = false;
%==========================================================================
% Plot parameters
%==========================================================================
sPlotParams.outputFolder                = 'figs';
sPlotParams.PlotEigenFuncsM             = 5;
sPlotParams.b_showEigenFigures          = false;
sPlotParams.b_showGraphMatricesFigures  = false;
sPlotParams.b_GSPBoxPlots               = true;
sPlotParams.b_plotSamplingPointsMarkers = false;
sPlotParams.b_plotHistograms            = false;

%% Generate graph
[G, ~, ~, sKernelParams] = GenerateGraph(graphName, nComponents, estDataDist, omega, true, M_G, N);
X = G.coords; % For convenience
if sPlotParams.b_plotHistograms
    graphNameSpaces = strrep(graphName,'_',' ');
    PlotHistogram(sPlotParams,X,graphNameSpaces,'$\hat{p}({\bf v})$ on $G$');
end

if sPlotParams.b_showGraphMatricesFigures
    PlotGraphMatrices(G, true);
end
if size(X,2) == 1
    sPlotParams.b_GSPBoxPlots = false;
end

[V, Lambda] = CalcNumericEigenvectors(M_G, sKernelParams, X);
if sPlotParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sPlotParams, estDataDist, X, [], 0, ...
        sPlotParams.PlotEigenFuncsM-1, V, Lambda, 'Numeric', [], G, ...
        '${\bf W}_G$ Eigenvectors')
end

phi = V(:,vInd);
PlotEigenfuncvecScatter(sPlotParams, estDataDist, X, [], v_m, v_m, V, Lambda, ...
    'Numeric', [], G, '${\bf W}_G$ Eigenvectors (Numeric)')

%% Nystrom
nysRatio = 0.8;
[VNys, LambdaNys] = CalcNystromEigenvectors(M_G, sKernelParams, X, nysRatio);
PlotEigenfuncvecScatter(sPlotParams, estDataDist, X, [], v_m, v_m, VNys, LambdaNys, ...
    'Nystrom', [], G, '${\bf W}_G$ Eigenvectors (Nystrom)')
%% G_tilde
%==========================================================================
% Transform vertices
%==========================================================================
if strcmp(R_TransformType, 'permutation')
    R = eye(G.N);
    r = randperm(G.N);
    R = R(r,:);
elseif strcmp(R_TransformType, 'eye')
    R = eye(G.N);
elseif strcmp(R_TransformType, 'randomMatrix')
    R = (1/sqrt(G.N))*randn(G.N,G.N);
end
X_tilde = R*X;

disp(['std(v_tilde)  = ' num2str(std(X_tilde))])
disp(['mean(v_tilde) = ' num2str(mean(X_tilde))])
disp(['std(R(:))     = ' num2str(std(R(:)))])
disp(['mean(R(:))    = ' num2str(mean(R(:)))])
disp(['1/sqrt(N)     = ' num2str(1/sqrt(N))])

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, ~, ~, sKernelTildeParams] = GenerateGraphTilde(X_tilde, nComponents,omega,M_Gtilde,true);
G_tilde_plt_title = sprintf('Gaussian %dD, %d components', size(X_tilde,2), nComponents);
if sPlotParams.b_plotHistograms
    PlotHistogram(sPlotParams,X_tilde, G_tilde_plt_title, '$\hat{p}(\tilde{{\bf v}})$ on $\tilde{G}$', true);
end

%==========================================================================
% Eigenfunctions
%==========================================================================
if strcmpi(phiTildeType, 'Analytic') || strcmpi(phiTildeType, 'Compare')
    [PhiAnalyticTilde, LambdaAnalyticTilde] = CalcAnalyticEigenfunctions(M_Gtilde, sKernelTildeParams, X_tilde, true);
    if sPlotParams.b_showEigenFigures
        firstEigenfunctionToPlot = 0;
        lastEigenfunctionToPlot = 4;
        PlotEigenfuncvecScatter(sPlotParams, estDataDist, X_tilde, [], firstEigenfunctionToPlot, lastEigenfunctionToPlot, ...
            PhiAnalyticTilde, LambdaAnalyticTilde, 'Analytic', [], G_tilde, ...
            ['${\bf W}_{\tilde{G}}$ kernel Eigenfunctions' newline '$\tilde{v}={\bf R} v$'])
        PlotInnerProductMatrix(sPlotParams, sKernelTildeParams.sDistParams, graphName, [], PhiAnalyticTilde, 'Analytic');
        PlotEigenvalues(sPlotParams, G_tilde_plt_title,'$\lambda_m$',Lambda, '$\tilde{\lambda}_m$', LambdaAnalyticTilde);
    end
    PhiTilde = PhiAnalyticTilde;
    LambdaTilde = LambdaAnalyticTilde;
    
elseif strcmpi(phiTildeType, 'Numeric')  || strcmpi(phiTildeType, 'Compare')
    [PhiNumericTilde, LambdaNumericTilde] = CalcNumericEigenvectors(M_Gtilde, sKernelTildeParams, X_tilde);
    if strcmpi(phiTildeType, 'Compare')
        PhiNumericTilde = FlipSign(PhiAnalyticTilde, PhiNumericTilde);
    end
    if sPlotParams.b_showEigenFigures
        firstEigenvectorToPlot = 0;
        lastEigenvectorToPlot = 4;
        PlotEigenfuncvecScatter(sPlotParams, estDataDist, X_tilde, [], firstEigenvectorToPlot, ...
            lastEigenvectorToPlot, PhiNumericTilde, LambdaNumericTilde, 'Numeric', [], G, ...
            '${\bf W}_G$ Eigenvectors')
        PlotEigenvalues(sPlotParams, G_tilde_plt_title,'$\lambda_m$',Lambda, '$\tilde{\lambda}_m$', LambdaNumericTilde);
    end
    PhiTilde = PhiNumericTilde;
    LambdaTilde = LambdaNumericTilde;
else
    error('invalid choise');
end
eigTh = 1e-10;
lastEigenvalueInd = find(LambdaTilde < eigTh, 1);
fprintf('*********************************************************\n');
fprintf('Last eigenvalue greater than %d is lambda #%d = %d\n', eigTh, lastEigenvalueInd-1, LambdaTilde(lastEigenvalueInd-1));
fprintf('*********************************************************\n');

if strcmpi(phiTildeType, 'Compare')
    figure;
    plot(vecnorm(PhiNumericTilde - PhiAnalyticTilde))
    title('$\| v_i - \phi_i \|$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
    keyboard;
end


%% Functional maps: C = L2.evecs'*L2.A'*P'*L1.evecs;
if strcmp(C_TransformType, 'PhiTilde^(-1)*R*Phi')
    if strcmpi(phiTildeType, 'Analytic')
        C = pinv(PhiTilde)*R*V;
    elseif strcmpi(phiTildeType, 'Numeric')
        C = PhiTilde'*R*V;
    end
    
    C_title = '\tilde{{\bf \Phi}}^\dagger {\bf R} {\bf \Phi}';
elseif strcmp(C_TransformType, 'PhiTilde^(-1)*Phi')
    if strcmpi(phiTildeType, 'Analytic')
        C = pinv(PhiTilde)*V;
    elseif strcmpi(phiTildeType, 'Numeric')
        C = PhiTilde'*V;
    end
    C_title = '\tilde{{\bf \Phi}}^\dagger {\bf \Phi}';
else
    error('select funcTransform');
end

PlotFmap(C,R,PhiTilde);


%% Transform to G tilde
Rphi = R*V(:,vInd); % R*PhiNumeric*alpha
phi_tilde = PhiTilde(:,vInd); % Phi_tilde(:,phiInd)
Rphi_title = ['${\bf R}\phi_{' num2str(vInd) '}$'];
phi_tilde_title = ['$\tilde{\phi_{' num2str(vInd) '}}$'];

if size(X_tilde,2) == 1
    figure;
    plot(X_tilde, R*V(:,vInd), 'o', 'DisplayName', Rphi_title)
    hold on
    plot(X_tilde, PhiTilde(:,vInd), '.', 'DisplayName', phi_tilde_title)
    legend('Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
else
    PlotGraphToGraphTransform(sPlotParams, G_tilde, '$\tilde{G}$',  G_tilde, '$\tilde{G}$',[],[],...
        Rphi, Rphi_title, phi_tilde, phi_tilde_title)    
end



%% Transform back to G
PhiNumeric_rec = R\PhiTilde*C; % R^(-1)*PhiAnalytic_tilde*C
% PhiNumeric_rec = Phi_tilde*C; % TODO: works only when size of C is NxN....
% phi_rec = Phi_tilde*C(:,phiInd); % Same as phi_rec = PhiNumeric_rec(:,phiInd);
% phi_rec = PhiNumeric*pinv(C)*alpha_tilde;
phi_rec = PhiNumeric_rec(:,vInd);
%==========================================================================
% Test phi_rec
%==========================================================================
fprintf('Phi%d\t\t Phi_rec%d\n', vInd, vInd);
disp([ V(1:5,vInd) phi_rec(1:5)])
%==========================================================================
% Plot transformation (without interpolation)
%==========================================================================
phi_title = ['$\phi_{' num2str(vInd) '}$ (Numeric)'];
f_tilde_title = ['$\tilde{f}$'];
phi_rec_title = ['$\phi_{{\bf rec}, ' num2str(vInd) '} $ (Numeric)'];
graphNameSpaces = strrep(graphName,'_',' ');
PlotGraphToGraphTransform(sPlotParams, G, ['$G$ (' graphNameSpaces ')'], G_tilde, '$\tilde{G}$', ...
    G, ['$G$ (' graphNameSpaces ')'], V(:,vInd), phi_title, f_tilde, f_tilde_title, phi_rec, phi_rec_title)
