%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
outputFolder = 'figs'; figSaveType = 'png';

%% Random data
signalType = 'Harmonics_1D'; % 'Harmonics_1D' / 'Audio_1D'
n = 2000;
if strcmp(signalType, 'Harmonics_1D')
    fs = 1000;
    Ts = 1/fs;
    xTrain = (0:n-1)'*Ts;
    
    f0 = 0.5; 
    f1 = 1;
    f2 = 3;
    signal = sin(2*pi*f0*xTrain) + sin(2*pi*f1*xTrain);% + sin(2*pi*f2*xTrain);
    xMax = max(xTrain);
    xMin = min(xTrain);
    
    fftLen = n; %2^ceil(log2(n));
    signalFft = fftshift(fft(signal, fftLen));
    fAxis = ((fs/fftLen)*(0:fftLen-1))' - fs/2; % [-fs/2, fs/2-fs/fftLen]
    
    omega = 0.3;
elseif strcmp(signalType, 'Audio_1D')
    [audio, fs] = audioread('Recording.m4a');
    n = 10000;
    audioLen = length(audio);
    signal = audio(round(audioLen*0.4):10:round(audioLen*0.6),1);
    soundsc(signal, fs);
    xTrain = (0:length(signal)-1)'/fs;

    xMax = max(xTrain);
    xMin = min(xTrain);
    
    fftLen = n; %2^ceil(log2(n));
    signalFft = fftshift(fft(signal, fftLen));
    fAxis = ((fs/fftLen)*(0:fftLen-1))' - fs/2; % [-fs/2, fs/2-fs/fftLen]
else
    error('invalid verticesPDF');
end

fig = figure('Name', 'Signal'); 
subplot(2,1,1)
    plot(xTrain,signal, '.')
    title('$s(t)$', 'interpreter', 'latex', 'FontSize', 16); 
    set(gca,'FontSize', 14);
subplot(2,1,2);
    stem(fAxis, abs(signalFft), 'filled');
    xlim([-1.5 1.5]*max([f0 f1 f2]));
    title('$S(f) = {\bf FFT}\{s(t)\}$', 'interpreter', 'latex', 'FontSize', 16); 
    set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig0_signal'), figSaveType);

fig = figure('Name', 'Histogram of X'); 
histogram(xTrain,100);
title('Histogram of ${\bf X}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig1_histogram_X'), figSaveType);
%% Generate graph
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
legend(strcat('$v_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig2_evecs_W'), figSaveType);

%% Project
alpha = pinv(V)*signal;
fig = figure('Name', 'Signal projection'); 
plot(xTrain, signal, 'o');
hold on
plot(xTrain, V*alpha, '.');
title('$s(t)$ vs. its projection', 'interpreter', 'latex', 'FontSize', 16);
legend('$s(t)$', '$\sum_{i=1}^M \alpha_i v_i(t)$', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
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
sigmaTilde = 0.5;
xTildeTrain = T(pCdf, true, muTilde, sigmaTilde, xTrain);

%% Demonstrate T (1/2)
xTestGrid = linspace(xMin,xMax,2000)';
polyCdf_xTestGrid = polyval(pCdf, xTestGrid);
b_saturate = false;
if b_saturate && (any(polyCdf_xTestGrid > 1) || any(polyCdf_xTestGrid < 0))
    warning('CDF must be in [0,1], saturating...');
    polyCdf_xTestGrid(polyCdf_xTestGrid > 1) = 1-eps; % saturate
    polyCdf_xTestGrid(polyCdf_xTestGrid < 0) = eps; % saturate
end

xTildeTestGrid = icdf('Normal',polyCdf_xTestGrid,muTilde,sigmaTilde);

% Generate the polynomial title
pCdfStr = [];
for p = 1:pCdfDegree+1
    if abs(pCdf(p)) < 1e-4
        continue;
    end
    if p > 1 && pCdf(p) > 0 && ~isempty(pCdfStr)
        pCdfStr = strcat(pCdfStr,'+');
    end
    if p == pCdfDegree+1
        pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree+1),'%.5f'));
    elseif p == pCdfDegree
        pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree),'%.5f'),'x');
    else
        pCdfStr = strcat(pCdfStr, num2str(pCdf(p),'%.5f'),'x^{',num2str(pCdfDegree-p+1),'}');
    end
    
end
fig = figure('Name', 'Demonstrate T (1/2)');
subplot(2,2,1)
    plot(xTrainGrid, estCdf_xTrainGrid,'.');
    hold on;
    plot(xTestGrid, polyCdf_xTestGrid, '.');
    ylim([0 1])
    title(['Est. vs. poly CDF (degree ' num2str(pCdfDegree) ')'], 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$\hat{F}_{X}(x)$', 'interpreter', 'latex', 'FontSize', 16);
    legend('${\bf eCDF}_{{\bf X}}(x_{{\bf train}})$', '$\hat{F}_{X}(x_{{\bf test}})$', ...
        'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(2,2,2)
    plot(xTildeTestGrid,polyCdf_xTestGrid,'.')
    title(['$\tilde{x} = T(x) = F_{\tilde{X}}^{-1}(\hat{F}_{X}(x))$' newline '(flipped axes)'], 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$\hat{F}_{X}(x)$', 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(2,2,3)
    histogram(xTestGrid ,100);
    title('Histogram of $x_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(2,2,4)
    histfit(xTildeTestGrid ,100);
    title('Histogram of $\tilde{x}_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
sgtitle(['$\hat{F}_{X}(x) = ' pCdfStr '$'], 'interpreter', 'latex', 'FontSize', 16);
saveas(fig,strcat(outputFolder, filesep, 'fig3_polyfit'), figSaveType);

%% Demonstrate T (2/2)
nTestPoints = 50;
xSmallTest = linspace(xMin+0.5,xMax-0.5,nTestPoints)';

xTestTilde = T(pCdf, true, muTilde, sigmaTilde, xSmallTest);
xSmallTestEst = invT(invpCdf, muTilde, sigmaTilde, xTestTilde);

cmap = xSmallTest;
fig = figure('Name', 'Demonstrate T (2/2)');
subplot(2,1,1)
    scatter(xSmallTest, zeros(1,nTestPoints), 100, cmap, 'o')
    hold on;
    scatter(xSmallTestEst, zeros(1,nTestPoints), 50, cmap, 'filled')
    colormap('jet')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    legend('$x_{{\bf test}}$', '$T^{-1}(T(x_{{\bf test}}))$','interpreter', 'latex', 'FontSize', 14);
    title('Original nodes','interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
subplot(2,1,2);
    scatter(xTestTilde, zeros(1,nTestPoints), 50, cmap, 'filled')
    colormap('jet')
    xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
    title('Transformed nodes','interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig4_x_to_xTilde'), figSaveType);

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
legend([strcat('$\tilde{v}_',string(vInd),'$') strcat('$\tilde{\phi}_',string(vInd),'$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
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
N = n*10;
xInt = linspace(xMin, xMax, N)';
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

%% Interpolated signal
figure('Name', 'Interpolated signal'); 
plot(xTrain, signal, 'o');
hold on
plot(xInt, VIntRenormed*alpha, '.');
title('Interpolated $s(t)$', 'interpreter', 'latex', 'FontSize', 16); 
legend('$s(t)$', '$\sum_{i=1}^M \alpha_i v_{{\bf int},i}(t)$', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
%% T(x)
function xTilde = T(pCdf, b_saturate, mu, sigma, x)
polyCdf = polyval(pCdf, x);
if (any(polyCdf > 1) || any(polyCdf < 0))
    warning('CDF must be in [0,1]...');
    if b_saturate
        polyCdf(polyCdf > 0.99) = 0.99; % saturate
        polyCdf(polyCdf < 0.01) = 0.01; % saturate
    end
end
xTilde = icdf('Normal',polyCdf,mu,sigma);
assert(~any(isnan(xTilde)),['xTilde contain NaNs since because of x = ', num2str(x(isnan(xTilde))')]);
end

%% invT(x)
function x = invT(invpCdf, mu, sigma, xTilde)
assert(~any(isnan(xTilde)),'xTilde contain NaNs...');
xTildeCdf = cdf('Normal', xTilde, mu, sigma);
x = polyval(invpCdf, xTildeCdf);
assert(~any(isnan(x)),'x contain NaNs...');
end