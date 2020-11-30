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
M_Gtilde = 10;

phiInd = 4;

funcTransform     = 'pinv(Btilde)RB'; % 'pinv(Btilde)RB' / 'pinv(Btilde)B'
verticesTransform = 'permutation'; % 'permutation' / 'randomMatrix' / 'eye'
%==========================================================================
% Plot parameters
%==========================================================================
sPlotParams.outputFolder                = 'figs';
sPlotParams.PlotEigenFuncsM             = 5;
sPlotParams.b_showEigenFigures          = false;
sPlotParams.b_showGraphMatricesFigures  = false;
sPlotParams.b_GSPBoxPlots               = true;
sPlotParams.b_plotSamplingPointsMarkers = true;
sPlotParams.b_plotHistograms            = true;

%% Generate graph
[G, ~, ~, sKernelParams] = GenerateGraph(graphName, nComponents, estDataDist, omega, true, M_G, N);
v = G.coords; % For convenience
if sPlotParams.b_plotHistograms
    graphNameSpaces = strrep(graphName,'_',' ');
    PlotHistogram(sPlotParams,v,graphNameSpaces,'$\hat{p}({\bf v})$ on $G$');
end

if sPlotParams.b_showGraphMatricesFigures
    PlotGraphMatrices(G, true);
end
if size(v,2) == 1
    sPlotParams.b_GSPBoxPlots = false;
end

[PhiNumeric, vLambdaNumeric] = CalcNumericEigenvectors(M_G, sKernelParams, v);
if sPlotParams.b_showEigenFigures
    PlotEigenfuncvecScatter(sPlotParams, estDataDist, v, [], 0, ...
        sPlotParams.PlotEigenFuncsM-1, PhiNumeric, vLambdaNumeric, 'Numeric', [], G, ...
        '${\bf W}_G$ Eigenvectors')
end

%% Generate graph signal
f = PhiNumeric(:,phiInd);
alpha = zeros(M_G,1);
alpha(phiInd) = 1;
%% G_tilde
%==========================================================================
% Transform
%==========================================================================
if strcmp(verticesTransform, 'permutation')
    R = eye(G.N);
    r = randperm(G.N);
    R = R(r,:);
elseif strcmp(verticesTransform, 'eye')
    R = eye(G.N);
elseif strcmp(verticesTransform, 'randomMatrix')
    R = (1/sqrt(G.N))*randn(G.N,G.N);
end
v_tilde = R*v;

disp(['std(v_tilde)  = ' num2str(std(v_tilde))])
disp(['mean(v_tilde) = ' num2str(mean(v_tilde))])
disp(['std(R(:))     = ' num2str(std(R(:)))])
disp(['mean(R(:))    = ' num2str(mean(R(:)))])
disp(['1/sqrt(N)     = ' num2str(1/sqrt(N))])

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, ~, ~, sKernelTildeParams] = GenerateGraphTilde(v_tilde, nComponents,omega,M_Gtilde,true);
if sPlotParams.b_plotHistograms
    G_tilde_plt_title = sprintf('Gaussian %dD, %d components', size(v_tilde,2), nComponents);
    PlotHistogram(sPlotParams,v_tilde, G_tilde_plt_title, '$\hat{p}(\tilde{{\bf v}})$ on $\tilde{G}$', true);
end

%==========================================================================
% Eigenfunctions
%==========================================================================
[PhiAnalytic_tilde, vLambdaAnalytic_tilde] = CalcAnalyticEigenfunctions(M_Gtilde, sKernelTildeParams, v_tilde, true);
if sPlotParams.b_showEigenFigures
    firstEigenfunctionToPlot = 0;
    lastEigenfunctionToPlot = 4;
    PlotEigenfuncvecScatter(sPlotParams, estDataDist, v_tilde, [], firstEigenfunctionToPlot, lastEigenfunctionToPlot, ...
        PhiAnalytic_tilde, vLambdaAnalytic_tilde, 'Analytic', [], G_tilde, ...
        ['${\bf W}_{\tilde{G}}$ kernel Eigenfunctions' newline '$\tilde{v}={\bf R} v$'])
    PlotInnerProductMatrix(sPlotParams, sKernelTildeParams.sDistParams, graphName, [], PhiAnalytic_tilde, 'Analytic');
end


[PhiNumeric_tilde, vLambdaNumeric_tilde] = CalcNumericEigenvectors(M_Gtilde, sKernelTildeParams, v_tilde);
PhiNumeric_tilde = FlipSign(PhiAnalytic_tilde, PhiNumeric_tilde);
if sPlotParams.b_showEigenFigures
    firstEigenvectorToPlot = 0;
    lastEigenvectorToPlot = 4;
    PlotEigenfuncvecScatter(sPlotParams, estDataDist, v_tilde, [], firstEigenvectorToPlot, ...
        lastEigenvectorToPlot, PhiNumeric_tilde, vLambdaNumeric_tilde, 'Numeric', [], G, ...
        '${\bf W}_G$ Eigenvectors')
end

figure; 
plot(vecnorm(PhiNumeric_tilde - PhiAnalytic_tilde))
title('$\| v_i - \phi_i \|$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

Phi_tilde = PhiAnalytic_tilde;
vLambda_tilde = vLambdaAnalytic_tilde;

lambda_title = '$\lambda_m$';
lambda_tilde_title = '$\tilde{\lambda}_m$';
PlotEigenvalues(sPlotParams, G_tilde_plt_title,lambda_title,vLambdaNumeric, lambda_tilde_title, vLambda_tilde);

%% Functional maps: C = L2.evecs'*L2.A'*P'*L1.evecs;
if strcmp(funcTransform, 'pinv(Btilde)RB')
    C = pinv(Phi_tilde)*R*PhiNumeric;
    C_title = '\tilde{{\bf \Phi}}^\dagger {\bf R} {\bf \Phi}';
elseif strcmp(funcTransform, 'pinv(Btilde)B')
    C = pinv(Phi_tilde)*PhiNumeric;
    C_title = '\tilde{{\bf \Phi}}^\dagger {\bf \Phi}';
else
    error('select funcTransform');
end

PlotFmap(C,R,Phi_tilde,phiInd);


%% Transform to G tilde
alpha_tilde = C*alpha; % C(:,phiInd)

figure;
subplot(1,2,1)
    scatter(1:M_G, alpha, 'bx', 'DisplayName', '$\alpha$')
    title('Coefficients on $G$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(1,2,2)
    scatter(1:M_Gtilde, alpha_tilde, 'ro', 'DisplayName', '$\tilde{\alpha}$')
    hold on
    scatter(1:M_Gtilde, C(:,phiInd), 'k+', 'DisplayName', ['${\bf C}(:,' num2str(phiInd) ')$'])
    title('Coefficients on $\tilde{G}$', 'Interpreter', 'latex', 'FontSize', 14)
    legend('Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
set(gcf,'Position', [400 400 1000 400])


% Cinv = pinv(C);
% f_tilde = PhiNumeric*Cinv(:,phiInd);
f_tilde = Phi_tilde*alpha_tilde;

Rphi = R*PhiNumeric(:,phiInd); % R*PhiNumeric*alpha
phi_tilde = Phi_tilde(:,phiInd); % Phi_tilde(:,phiInd)
Rphi_title = ['${\bf R}\phi_{' num2str(phiInd) '}$'];
phi_tilde_title = ['$\tilde{\phi_{' num2str(phiInd) '}}$'];

if size(v_tilde,2) == 1
    figure;
    plot(v_tilde, R*PhiNumeric(:,phiInd), 'o', 'DisplayName', Rphi_title)
    hold on
    plot(v_tilde, Phi_tilde(:,phiInd), '.', 'DisplayName', phi_tilde_title)
    legend('Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
else
    PlotGraphToGraphTransform(sPlotParams, G_tilde, '$\tilde{G}$',  G_tilde, '$\tilde{G}$',[],[],...
        Rphi, Rphi_title, phi_tilde, phi_tilde_title)    
end



%% Transform back to G
PhiNumeric_rec = R\Phi_tilde*C; % R^(-1)*PhiAnalytic_tilde*C
% PhiNumeric_rec = Phi_tilde*C; % TODO: works only when size of C is NxN....
% phi_rec = Phi_tilde*C(:,phiInd); % Same as phi_rec = PhiNumeric_rec(:,phiInd);
% phi_rec = PhiNumeric*pinv(C)*alpha_tilde; %okay
phi_rec = PhiNumeric_rec(:,phiInd);
%==========================================================================
% Test phi_rec
%==========================================================================
fprintf('Phi%d\t\t Phi_rec%d\n', phiInd, phiInd);
disp([ PhiNumeric(1:5,phiInd) phi_rec(1:5)])

% recvecnorm = vecnorm( PhiNumeric_rec - PhiNumeric, 2 );
% fprintf('First 5 values of vecnorm( PhiNumeric_rec - PhiNumeric, 2 )\n');
% disp(recvecnorm(1:5))

%==========================================================================
% Plot transformation (without interpolation)
%==========================================================================
phi_title = ['$\phi_{' num2str(phiInd) '}$ (Numeric)'];
% f_tilde_title = ['$\tilde{f} = \tilde{{\bf \Phi}} {\bf C} {\bf \alpha}$' newline ...
%     '$\big({\bf C} = ' C_title '\big)$'];
f_tilde_title = ['$\tilde{f}$'];
phi_rec_title = ['$\phi_{{\bf rec}, ' num2str(phiInd) '} $ (Numeric)'];
graphNameSpaces = strrep(graphName,'_',' ');
PlotGraphToGraphTransform(sPlotParams, G, ['$G$ (' graphNameSpaces ')'], G_tilde, '$\tilde{G}$', ...
    G, ['$G$ (' graphNameSpaces ')'], PhiNumeric(:,phiInd), phi_title, f_tilde, f_tilde_title, phi_rec, phi_rec_title)



%% Interpolate on Gtilde
if size(v,2) == 1
    % TODO: should think on how to interpolate the nodes on Gtilde...
    v_tilde_interp = interp(v_tilde,2);
    
%     figure; plot(v_tilde, ones(N,1), 'o'); hold on; plot(v_tilde_interp, ones(2*N,1), 'x')
    
    [PhiNumeric_tilde_interp, ~] = CalcNumericEigenvectors(M_Gtilde, sKernelTildeParams, v_tilde_interp);
    phi_interp_title = ['$\phi_{{\bf int}, ' num2str(phiInd) '}$ (Numeric)'];
    G_tilde_interp = G_tilde;
    G_tilde_interp.coords = v_tilde_interp;

%     f_tilde_interp = PhiNumeric_tilde_interp*alpha_tilde;
%     f_tilde_interp_title = ['$\tilde{f} = \tilde{{\bf \Phi}}_{{\bf int}} {\bf C} \tilde{{\bf \alpha}}$'];

    %==========================================================================
    % Transform back to G
    %==========================================================================
    v_interp = interp(v,2);
    G_interp = G;
    G_interp.coords = v_interp;

    phi_rec_interp = PhiNumeric_tilde_interp*C(:,phiInd);
    phi_rec_interp_title = ['$\phi_{{\bf rec}, ' num2str(phiInd) '} $ (Numeric)'];
    
    
    graphNameSpaces = strrep(graphName,'_',' ');
    PlotGraphToGraphTransform(sPlotParams, ...
        G, ['$G$ (' graphNameSpaces ')'], ...
        G_interp, ['$G_{{\bf int}}$ (' graphNameSpaces ')'], ...
        [],[],...
        PhiNumeric(:,phiInd), phi_title, ...
        phi_rec_interp, phi_rec_interp_title)

    
end

%% Extra debug plots


%==========================================================================
% Test alpha_rec
%==========================================================================
alpha_rec = pinv(C)*alpha_tilde;

figure; 
title('$\alpha$ vs. $\alpha_{{\bf rec}$', 'Interpreter', 'latex', 'FontSize', 14)
scatter(1:M_G, alpha, 'bx', 'DisplayName', '$\alpha$')
hold on
scatter(1:M_G, alpha_rec, 'ro', 'DisplayName', '$\alpha_{{\bf rec}}$')
legend('Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

fprintf('alpha_rec\t\talpha\n');
disp([alpha_rec(1:5) alpha(1:5)]);