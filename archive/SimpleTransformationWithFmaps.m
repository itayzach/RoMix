%% Restart
clc; clear; rng('default')
close all; 
%% Parameters
graphName =  'Uniform_2D'; % Uniform_1D
nNodes = 1000; % number of nodes in G

%==========================================================================
% kernel parameters
%==========================================================================
omega       = 0.3;
estDataDist = 'Gaussian';
nComponents = 1;

%==========================================================================
% basis parameters
%==========================================================================
M_G = 30;
M_Gtilde = 30;

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
sPlotParams.b_showVerticesTransform     = false;
sPlotParams.b_GSPBoxPlots               = true;
sPlotParams.b_plotTransformation        = false;
sPlotParams.b_calcalphaInterpolated     = false;
sPlotParams.b_interpolateOnGWithLS      = true;
sPlotParams.b_plotSamplingPointsMarkers = true;

%% Generate graph
[G, ~, ~, sKernelParams] = GenerateGraph(graphName, nComponents, estDataDist, omega, true, M_G, nNodes);
v = G.coords; % For convenience

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
alpha = 10*randn(M_G, 1);
% alpha = zeros(M_G, 1);
% alpha(phiInd) = 1;
f = PhiNumeric*alpha;

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

%==========================================================================
% Generate G_tilde
%==========================================================================
[G_tilde, ~, ~, sKernelTildeParams] = GenerateGraphTilde(v_tilde, nComponents,omega,M_Gtilde,true);

%==========================================================================
% Eigenfunctions
%==========================================================================
[PhiAnalytic_tilde, vLambdaAnalytic_tilde] = CalcAnalyticEigenfunctions(M_Gtilde, sKernelTildeParams, v_tilde, true);
if sPlotParams.b_showEigenFigures
    firstEigenfunctionToPlot = 144;
    lastEigenfunctionToPlot = 149;
    PlotEigenfuncvecScatter(sPlotParams, estDataDist, v_tilde, [], firstEigenfunctionToPlot, lastEigenfunctionToPlot, ...
        PhiAnalytic_tilde, vLambdaAnalytic_tilde, 'Analytic', [], G_tilde, ...
        ['${\bf W}_{\tilde{G}}$ kernel Eigenfunctions' newline '$\tilde{v}={\bf R} v$'])
    PlotInnerProductMatrix(sPlotParams, sKernelTildeParams.sDistParams, graphName, [], PhiAnalytic_tilde, 'Analytic');
end

% [PhiNumeric_tilde, vLambdaNumeric_tilde] = CalcNumericEigenvectors(M_Gtilde, sKernelTildeParams, v_tilde);
% % PhiNumeric_tilde = FlipSign(PhiAnalytic_tilde, PhiNumeric_tilde);
% if sPlotParams.b_showEigenFigures
%     firstEigenvectorToPlot = 144;
%     lastEigenvectorToPlot = 149;
%     PlotEigenfuncvecScatter(sPlotParams, estDataDist, v_tilde, [], firstEigenvectorToPlot, ...
%         lastEigenvectorToPlot, PhiNumeric_tilde, vLambdaNumeric_tilde, 'Numeric', [], G, ...
%         '${\bf W}_G$ Eigenvectors')
% end

% figure;
% idx = 10;
% plot(vecnorm(PhiNumeric_tilde(:,1:idx)-PhiAnalytic_tilde(:,1:idx)))

% Phi_tilde = PhiNumeric_tilde;
Phi_tilde = PhiAnalytic_tilde;

%% Functional maps: C = L2.evecs'*L2.A'*P'*L1.evecs;
if strcmp(funcTransform, 'pinv(Btilde)RB')
    C = pinv(Phi_tilde)*R*PhiNumeric;
elseif strcmp(funcTransform, 'pinv(Btilde)B')
    C = pinv(Phi_tilde)*PhiNumeric;
    T = Phi_tilde*pinv(Phi_tilde)*PhiNumeric*pinv(PhiNumeric);
else
    error('select funcTransform');
end
% T = Phi_tilde*C*pinv(PhiNumeric);
% 
alpha_tilde = C*alpha;
% f_tilde = T*f;
% f_tilde_title = '$\tilde{f} = {\bf T} f$';
% 
f_tilde = Phi_tilde*alpha_tilde;
f_tilde_title = '$\tilde{f}(\tilde{v}) = \tilde{{\bf \Phi}} \tilde{\alpha}$';

% 
% f_tilde3 = Phi_tilde*C;
% f_tilde3_title = '$\tilde{f} = \tilde{{\bf \Phi}} \tilde{\bf C}$';
%% Test
%==========================================================================
% Show C, pinv(C) and R
%==========================================================================
figure; 
subplot(1,3,1)
    imagesc(R); colorbar;
    title('${\bf R}$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(1,3,2)
    imagesc(C); colorbar;
    title('${\bf C}$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(1,3,3)    
    imagesc(pinv(C)); colorbar;
    title('${\bf C}^\dagger$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
set(gcf,'Position', [400 400 1500 400])


% %==========================================================================
% % Test identity
% %==========================================================================
% PhiTilde_PhiTildeInv = Phi_tilde*pinv(Phi_tilde);
% PhiTildeInv_PhiTilde = pinv(Phi_tilde)*Phi_tilde;
% 
% figure; 
% subplot(1,2,1)
%     imagesc(PhiTilde_PhiTildeInv); colorbar;
%     title('$\tilde{\Phi} \Phi^\dagger$', 'Interpreter', 'latex', 'FontSize', 14)
%     set(gca,'FontSize', 14);
% subplot(1,2,2)
%     imagesc(PhiTildeInv_PhiTilde); colorbar;
%     title('$\Phi^\dagger \tilde{\Phi}$', 'Interpreter', 'latex', 'FontSize', 14)
%     set(gca,'FontSize', 14);
% set(gcf,'Position', [400 400 1000 400])


if strcmp(funcTransform, 'pinv(Btilde)RB')
    f_rec = R\Phi_tilde*C*alpha;
    f_rec_title = '$f_{{\bf rec}}(v) = {\bf R}^{-1} \tilde{{\bf \Phi}} {\bf C} \alpha$';
elseif strcmp(funcTransform, 'pinv(Btilde)B')
    f_rec = Phi_tilde*C*alpha;
    f_rec_title = '$f_{{\bf rec}}(v) = \tilde{{\bf \Phi}} {\bf C} \alpha$';
else
    error('select funcTransform');
end



f_title = '$f(v) = {\bf \Phi} \alpha$';
PlotGraphToGraphTransform(sPlotParams, G, '$G$', G_tilde, '$\tilde{G}$', ...
    G, '$G$', f, f_title, f_tilde, f_tilde_title, f_rec, f_rec_title)


% %==========================================================================
% % Test alpha_rec
% %==========================================================================
% alpha_rec = pinv(C)*alpha_tilde;
% 
% figure; 
% title('$\alpha$ vs. $\alpha_{{\bf rec}$', 'Interpreter', 'latex', 'FontSize', 14)
% scatter(1:M_G, alpha, 'bx', 'DisplayName', '$\alpha$')
% hold on
% scatter(1:M_G, alpha_rec, 'ro', 'DisplayName', '$\alpha_{{\bf rec}}$')
% scatter(1:M_Gtilde, C(:,phiInd), 'k+', 'DisplayName', ['${\bf C}(:,' num2str(phiInd) ')$'])
% legend('Interpreter', 'latex', 'FontSize', 14)
% set(gca,'FontSize', 14);
% 
% fprintf('alpha_rec\t\talpha\n');
% disp([alpha_rec(1:5) alpha(1:5)]);

% %==========================================================================
% % Test phi_rec
% %==========================================================================
% % figure; plot(v, [Phi_tilde*C(:,phiInd), PhiNumeric(:,phiInd)],'.')
% 
% phi = PhiNumeric(:,phiInd);
% phi_title = ['$\phi_{' num2str(phiInd) '}$ (Numeric)'];
% 
% %     phi_tilde = Phi_tilde(:,phiInd);
% %     phi_tilde_title = ['$\tilde{\phi_{' num2str(phiInd) '}}$ (Analytic)'];
% 
% 
% % PhiNumeric_rec = R\Phi_tilde*C; % R^(-1)*PhiAnalytic_tilde*C
% PhiNumeric_rec = Phi_tilde*C; % R^(-1)*PhiAnalytic_tilde*C
% phi_rec = PhiNumeric_rec(:,phiInd);
% phi_rec_title = ['$\phi_{{\bf rec}, ' num2str(phiInd) '} $ (Numeric)'];
% 
% fprintf('Phi%d\t\t Phi_rec%d\n', phiInd, phiInd);
% disp([ phi(1:5) phi_rec(1:5)])
% 
% vecnorm( PhiNumeric_rec - PhiNumeric, 2 )
% 
% PlotGraphToGraphTransform(sPlotParams, G, '$G$', G_tilde, '$\tilde{G}$', ...
%     G, '$G$', phi, phi_title, f_tilde, f_tilde_title, phi_rec, phi_rec_title)
% 
% PlotGraphToGraphTransform(sPlotParams, G, '$G$', G_tilde, '$\tilde{G}$', ...
%     G, '$G$', phi, phi_title, f_tilde2, f_tilde2_title, phi_rec, phi_rec_title)
% 
% PlotGraphToGraphTransform(sPlotParams, G, '$G$', G_tilde, '$\tilde{G}$', ...
%     G, '$G$', phi, phi_title, f_tilde3, f_tilde3_title, phi_rec, phi_rec_title)


% %==========================================================================
% % Test f_rec
% %==========================================================================
% f_rec = pinv(T)*f_tilde;
% fprintf('f_rec\t\tf\n');
% disp([f_rec(1:5) f(1:5)]);
% f_title = '$f(v) = {\bf V} \alpha$';
% f_tilde_title = '$\tilde{f}(\tilde{v}) = \tilde{{\bf \Phi}} \tilde{\alpha}$';
% f_rec_title = '$f_{{\bf rec}}(v) = \tilde{{\bf \Phi}} \alpha_{{\bf rec}}$';
% 
% PlotGraphToGraphTransform(sPlotParams, G, '$G$', G_tilde, '$\tilde{G}$', ...
%     G, '$G$', f, f_title, f_tilde, f_tilde_title, f_rec, f_rec_title)

%% Interpolate on Gtilde
% if size(v,2) == 1
%     % TODO: should think on how to interpolate the nodes on Gtilde...
%     v_tilde_interp = interp(v_tilde,2);
% 
%     [BtildeInterp, ~] = CalcAnalyticEigenfunctions(M_Gtilde, sKernelTildeParams, ...
%             v_tilde_interp, true);
% 
%     G_tilde_interp = G_tilde;
%     G_tilde_interp.coords = v_tilde_interp;
% 
%     f_tilde_interp = BtildeInterp*alpha_tilde;
% 
%     %==========================================================================
%     % Transform back to G
%     %==========================================================================
%     %TODO: How to extend f_interp to 1000x1??? it's 500x1....
%     Tinterp = BtildeInterp*C*pinv(PhiNumeric);
%     f_interp = pinv(Tinterp)*f_tilde_interp;
%     f_interp_title = '$f_{{\bf int}}(v) = {\bf T}^\dagger \tilde{f}_{{\bf int}}(\tilde{v})$';
%     f_title = '$f(v) = {\bf V} \alpha$';
%     f_tilde_interp_title = '$\tilde{f}_{{\bf int}}(\tilde{v}) = \tilde{{\bf \Phi}} \tilde{c}$';
% 
%     PlotGraphToGraphTransform(sPlotParams, G, '$G$', G_tilde_interp, '$\tilde{G}$', ...
%         G, '$G$', f, f_title, f_tilde_interp, f_tilde_interp_title, f_interp, f_interp_title)
% end