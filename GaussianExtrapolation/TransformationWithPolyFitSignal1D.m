%% Restart
clc; clear; close all; rng('default');
set(0,'DefaultFigureWindowStyle','docked')
outputFolder = 'figs'; figSaveType = 'png';

%% Random data
signalType = 'Synthetic_1D'; % 'Synthetic_1D' / 'Audio_1D'
n = 1000;

fs = 48e3;
Ts = 1/fs;
tTrain = Ts*(0:n-1)';
tMax = max(tTrain);
tMin = min(tTrain);
if strcmp(signalType, 'Synthetic_1D')
    f1 = 1*0.5/(n*Ts); % m*(0.5fs/n)
    f2 = 3*0.5/(n*Ts);
    f3 = 5*0.5/(n*Ts);

    signalTrain = sin(2*pi*f1*tTrain) + sin(2*pi*f2*tTrain) + sin(2*pi*f3*tTrain);
elseif strcmp(signalType, 'Audio_1D')
    [signalTrain, fs_audioread] = audioread('data/signal.m4a');
    assert(fs == fs_audioread);
    signalTrain = signalTrain(:,1); % mono-channel
    signalTrain = signalTrain(abs(signalTrain) > 0); % take only samples with energy
    signalTrain = signalTrain(1:n); % take only first n samples (n*Ts=20ms for fs=48KHz, n = 1000)
    % filter
    windowSize = 100; 
    b = ones(1,windowSize);
    a = 1;
    signalTrain = (1/windowSize)*filter(b,a,signalTrain);
%     soundsc(signal, fs);
    
else
    error('invalid verticesPDF');
end

fftLen = n;
signalFft = fftshift(fft(signalTrain, fftLen));
freqAxis = ((fs/fftLen)*(0:fftLen-1))' - fs/2; % [-fs/2, fs/2-fs/fftLen]

fig = figure('Name', 'Signal'); 
subplot(2,1,1)
    plot(tTrain*1e3,signalTrain, 'LineWidth', 2)
    title('$s(t)$', 'interpreter', 'latex', 'FontSize', 16); 
    xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14); 
    legend('$s_{{\bf train}}(t)$', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(2,1,2);
    plot(freqAxis/1e3, 20*log10(abs(signalFft)));
    xlim([-0.5 0.5]*fs/1e3);
    title('$S(f) = {\bf FFT}\{s(t)\}$', 'interpreter', 'latex', 'FontSize', 16); 
    xlabel('$f$ [KHz]', 'interpreter', 'latex', 'FontSize', 14); 
    ylabel('$S(f)$ [dB]', 'interpreter', 'latex', 'FontSize', 14); 
    set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig0_signal'), figSaveType);

fig = figure('Name', 'Histogram of tTrain'); 
histogram(tTrain*1e3,100);
title('Histogram of $t_{{\bf train}}$', 'interpreter', 'latex', 'FontSize', 16); 
xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig1_histogram_X'), figSaveType);
%% Generate graph
omega = Ts;
dist = pdist2(tTrain, tTrain);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
M = 20;
% d = sum(W,1);
% Dn = diag(d.^(-0.5));
% L = eye(n) - Dn * W * Dn;
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:5;

fig = figure('Name', 'Eigenvectors of W');
plot(tTrain*1e3, V(:,vInd),'.');
title('Eigenvectors of ${\bf W}$ ', 'interpreter', 'latex', 'FontSize', 16);
xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
legend(strcat('$v_{',string(vInd),'}$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig2_evecs_W'), figSaveType);


VFft = fftshift(fft(V, fftLen));   
figure('Name', 'FFT of evecs of W');
plot(freqAxis/1e3, 20*log10(abs(VFft(:,vInd))),'LineWidth',2);
title('FFT of eigenvectors', 'interpreter', 'latex', 'FontSize', 16); 
legend(strcat('$V_{',string(vInd),'}(f)$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
xlabel('$f$ [KHz]', 'interpreter', 'latex', 'FontSize', 14); 
ylabel('$S(f)$ [dB]', 'interpreter', 'latex', 'FontSize', 14); 
set(gca,'FontSize', 14);
%% Project
alpha = pinv(V)*signalTrain;
projSignalFft = fftshift(fft(V*alpha, fftLen));   
figure('Name', 'Signal projection'); 
subplot(2,1,1)
    plot(tTrain*1e3, signalTrain, 'o');
    hold on
    plot(tTrain*1e3, V*alpha, '.');
    title(['$s(t)$ vs. its projection with $M = ' num2str(M) '$'], 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
    legend('$s(t)$', '$\sum_{i=1}^M \alpha_i v_i(t)$', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(2,1,2);
    plot(freqAxis/1e3, 20*log10(abs(signalFft)));
	xlim([-0.5 0.5]*fs/1e3);
    hold on;
    plot(freqAxis/1e3, 20*log10(abs(projSignalFft)));
    title('$S(f) = {\bf FFT}\{s(t)\}$', 'interpreter', 'latex', 'FontSize', 16); 
    legend('$S(f)$', '$\sum_{i=1}^M \alpha_i V_i(f)$', 'interpreter', 'latex', 'FontSize', 14)
    xlabel('$f$ [KHz]', 'interpreter', 'latex', 'FontSize', 14); 
    ylabel('$S(f)$ [dB]', 'interpreter', 'latex', 'FontSize', 14); 
    set(gca,'FontSize', 14);
%% Learn CDF(t) from xTrain by fitting ecdf to a polynomial
[estCdf_tTrainGrid, tTrainGrid] = ecdf(tTrain);
tTrainGrid = tTrainGrid(1:end-1);
estCdf_tTrainGrid = estCdf_tTrainGrid(1:end-1);

pCdfDegree = 2;
invpCdfDegree = 2;
pCdf = polyfit(tTrainGrid, estCdf_tTrainGrid, pCdfDegree); % to have analytic expression for the cdf
invpCdf = polyfit(estCdf_tTrainGrid, tTrainGrid, invpCdfDegree); % to have analytic expression for the cdf

%% Transform to Gtilde
muTilde = 0;
sigmaTilde = 3*Ts;%(xMax-xMin)/sqrt(12);
tTildeTrain = T(pCdf, true, muTilde, sigmaTilde, tTrain);

%% Demonstrate T (1/2)
tTestGrid = linspace(tMin,tMax,2000)';
polyCdf_tTestGrid = polyval(pCdf, tTestGrid);
b_saturate = false;
if b_saturate && (any(polyCdf_tTestGrid > 1) || any(polyCdf_tTestGrid < 0))
    warning('CDF must be in [0,1], saturating...');
    polyCdf_tTestGrid(polyCdf_tTestGrid > 1) = 1-eps; % saturate
    polyCdf_tTestGrid(polyCdf_tTestGrid < 0) = eps; % saturate
end

tTildeTestGrid = icdf('Normal',polyCdf_tTestGrid,muTilde,sigmaTilde);

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
        pCdfStr = strcat(pCdfStr, num2str(pCdf(pCdfDegree),'%.5f'),'t');
    else
        pCdfStr = strcat(pCdfStr, num2str(pCdf(p),'%.5f'),'t^{',num2str(pCdfDegree-p+1),'}');
    end
    
end
fig = figure('Name', 'Demonstrate T (1/2)');
subplot(2,2,1)
    plot(tTrainGrid*1e3, estCdf_tTrainGrid,'.');
    hold on;
    plot(tTestGrid*1e3, polyCdf_tTestGrid, '.');
    ylim([0 1])
    title(['Est. vs. poly CDF (degree ' num2str(pCdfDegree) ')'], 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
    ylabel('$\hat{F}_{X}(t)$', 'interpreter', 'latex', 'FontSize', 16);
    legend('${\bf eCDF}_{{\bf X}}(t_{{\bf train}})$', '$\hat{F}_{X}(t_{{\bf test}})$', ...
        'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(2,2,2)
    plot(tTildeTestGrid*1e3,polyCdf_tTestGrid,'.')
    title(['$\tilde{t} = T(t) = F_{\tilde{X}}^{-1}(\hat{F}_{X}(t))$' newline '(flipped axes)'], 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$\hat{F}_{X}(t)$', 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$\tilde{t}$ [ms]', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(2,2,3)
    histogram(tTestGrid*1e3 ,100);
    title('Histogram of $t_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
subplot(2,2,4)
    histfit(tTildeTestGrid*1e3 ,100);
    title('Histogram of $\tilde{t}_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
    xlabel('$\tilde{t}$ [ms]', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
sgtitle(['$\hat{F}_{X}(t) = ' pCdfStr '$'], 'interpreter', 'latex', 'FontSize', 16);
saveas(fig,strcat(outputFolder, filesep, 'fig3_polyfit'), figSaveType);

%% Demonstrate T (2/2)
nTestPoints = 80;
tSmallTest = linspace(tTrain(2),tTrain(end-1),nTestPoints)';

tTestTilde = T(pCdf, true, muTilde, sigmaTilde, tSmallTest);
tSmallTestEst = invT(invpCdf, muTilde, sigmaTilde, tTestTilde);

cmap = tSmallTest;
fig = figure('Name', 'Demonstrate T (2/2)');
subplot(2,1,1)
    scatter(tSmallTest*1e3, zeros(1,nTestPoints), 100, cmap, 'o')
    hold on;
    scatter(tSmallTestEst*1e3, zeros(1,nTestPoints), 50, cmap, 'filled')
    colormap('jet')
    xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 16);
    legend('$x_{{\bf test}}$', '$T^{-1}(T(x_{{\bf test}}))$','interpreter', 'latex', 'FontSize', 14);
    title('Original nodes','interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
subplot(2,1,2);
    scatter(tTestTilde*1e3, zeros(1,nTestPoints), 50, cmap, 'filled')
    colormap('jet')
    xlabel('$\tilde{t}$ [ms]', 'interpreter', 'latex', 'FontSize', 16);
    title('Transformed nodes','interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig4_x_to_xTilde'), figSaveType);

%% Build G tilde
omegaTilde = omega;
distTilde = pdist2(tTildeTrain, tTildeTrain);
WTilde = exp(-distTilde.^2/(2*omegaTilde^2));

%% Calculate analytic eigenfunctions of W_tilde, and numeric eigenvectors for comparison
MTilde = 20;
[PhiTilde, ~] = SimpleCalcAnalyticEigenfunctions(tTildeTrain, omegaTilde, sigmaTilde, muTilde, MTilde);
[VTilde, LambdaNumericTilde] = eigs(WTilde, MTilde);
VTilde = FlipSign(PhiTilde, VTilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

fig = figure('Name', 'evecs vs. efuncs of WTilde');
plot(tTildeTrain*1e3, VTilde(:,vInd),'o');
hold on
plot(tTildeTrain*1e3, PhiTilde(:,vInd),'.');
xlabel('$\tilde{t}$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
legend([strcat('$\tilde{v}_{',string(vInd),'}$') strcat('$\tilde{\phi}_{',string(vInd),'}$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
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
tEstTrain = invT(invpCdf, muTilde, sigmaTilde, tTildeTrain);
VRec = PhiTilde*C;

fig = figure('Name', 'V in terms of PhiTilde');
plot(tTrain*1e3, V(:,vInd),'o');
hold on
plot(tEstTrain*1e3, VRec(:,vInd),'.');
xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
legend([strcat('$v_{',string(vInd),'}$') strcat('$v_{{\bf rec},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['${\bf V}$ in terms of ${\bf \tilde{\Phi}}$' newline '${\bf V}_{{\bf rec}} = {\bf \tilde{\Phi}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig8_VRec'), figSaveType);

%% Interpolate
interpRatio = 10;
N = interpRatio*n;
tInt = linspace(tTrain(10), tTrain(end-10), N)'; % cut edges
tTildeInt = T(pCdf, true, muTilde, sigmaTilde, tInt);
[PhiTildeInt, ~] = SimpleCalcAnalyticEigenfunctions(tTildeInt, omegaTilde, sigmaTilde, muTilde, MTilde);

fig = figure('Name', 'Eigenfunctions of WTilde on the entire axis');
plot(tTildeInt*1e3, PhiTildeInt(:,vInd),'LineWidth',2);
legend(strcat('$\tilde{\phi}_{',string(vInd),'}$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Eigenfunctions of $\tilde{{\bf W}}$ on the entire axis' newline ...
    '$\tilde{\omega}$ = ' num2str(omegaTilde*1e3, '%.2f') ...
    '[ms]; $\tilde{\sigma}$ = ' num2str(sigmaTilde*1e3, '%.2f') ...
    '[ms]; $\tilde{\mu}$ = ' num2str(muTilde*1e3, '%.2f') '[ms]'], 'interpreter', 'latex', 'FontSize', 16); 
xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig9_efuncs_WTilde_int'), figSaveType);

tIntInvT = invT(invpCdf, muTilde, sigmaTilde, tTildeInt);
VInt = PhiTildeInt*C;

VIntRenormed = sqrt(interpRatio)*VInt;

fig = figure('Name', 'Interpolated evecs of W');
plot(tTrain*1e3, V(:,vInd),'o');
hold on
plot(tIntInvT*1e3, VIntRenormed(:,vInd),'.');
xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);
legend([strcat('$v_{',string(vInd),'}$') strcat('$v_{{\bf int},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Interpolated eigenvectors of ${\bf W}$' newline '${\bf V}_{{\bf int}} = \sqrt{\frac{N}{n}}{\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
saveas(fig,strcat(outputFolder, filesep, 'fig10_VInt'), figSaveType);

%% Interpolated signal
intSignal = VIntRenormed*alpha;
intFftLen = N;
intSignalFft = (fftLen/intFftLen)*fftshift(fft(intSignal, intFftLen));
intFreqAxis = ((interpRatio*fs/intFftLen)*(0:intFftLen-1))' - interpRatio*fs/2; % [-fs/2, fs/2-fs/fftLen]

fig = figure('Name', 'Interpolated signal'); 
subplot(2,1,1)
    plot(tTrain*1e3, signalTrain, 'o');
    hold on
    plot(tIntInvT*1e3, intSignal, '.')
    title('Interpolated $s(t)$', 'interpreter', 'latex', 'FontSize', 16); 
    xlabel('$t$ [ms]', 'interpreter', 'latex', 'FontSize', 14);     
    legend('$s(t)$', '$\sum_{i=1}^M \alpha_i v_{{\bf int},i}(t)$', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(2,1,2);
    plot(freqAxis/1e3, 20*log10(abs(signalFft)));
    hold on
    plot(intFreqAxis/1e3, 20*log10(abs(intSignalFft)));
    xlim([-0.5 0.5]*fs/1e3);
    title('$S(f) = {\bf FFT}\{s(t)\}$', 'interpreter', 'latex', 'FontSize', 16); 
    legend('$S(f)$', '$\sum_{i=1}^M \alpha_i V_{{\bf int},i}(f)$', 'interpreter', 'latex', 'FontSize', 14)
    xlabel('$f$ [KHz]', 'interpreter', 'latex', 'FontSize', 14); 
    ylabel('$S(f)$ [dB]', 'interpreter', 'latex', 'FontSize', 14); 
    set(gca,'FontSize', 14);
