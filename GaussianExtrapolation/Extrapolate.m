function [] = Extrapolate(sParams, sSimParams)

dx = 0.01;
x = (-5:dx:5-dx)';
% x = (50:dx:60-dx)';
x = repmat(x, 1, sParams.dim);
assert(size(x, 2) == sParams.dim)
N = length(x);



A1 = 5;
B1 = 0.001;
C1 = 0.001;
A2 = 10;
B2 = 7;
C2 = 6;
% mF      = [ A*exp(-x).*(sin(2.5*x) + sin(2*pi*x))   A*exp(-2*x).*sin(5*x) ];
[phi1, ~] = phi(sParams.a, sParams.b, 1, x);
[phi2, ~] = phi(sParams.a, sParams.b, 6, x);
[phi3, ~] = phi(sParams.a, sParams.b, 10, x);
mF      = [ A1*phi1 + B1*phi2 + C1*phi3   ...
            A2*exp(-0.2*x.^2).*sin(pi*x) + B2*exp(-0.5*x.^2).*sin(2*pi*x) + B2*exp(-0.3*x.^2).*sin(1*pi*x) ];
mF_awgn = [sqrt(sSimParams.noiseVar1)*randn(N,1) sqrt(sSimParams.noiseVar2)*randn(N,1)];
% cFstr   = {'10e^{-x}\big(\sin(2.5x) + \sin(2\pi x)\big)' '10e^{-2x}\sin(5x)'};
nFuncs  = size(mF, 2);

% x_eval = ([-1:0.1:1])';
% vF_eval = [10*exp(-2*x).*sin(5*x)];

cFigs = cell(1, nFuncs);

for i = 1:nFuncs
    vFi = mF(:, i);
    vFi_awgn = mF_awgn(:, i);
    SNR = snr(vFi, vFi_awgn);
    vGi = vFi + vFi_awgn;
    mPhi = zeros(N, sParams.M);
    for m = 0:sParams.M-1 
        [vPhi_m_x, lambda_m] = phi(sParams.a, sParams.b, m, x);
        mPhi(:, m+1) = vPhi_m_x;
    end
    assert(sParams.M <= sParams.R, 'You cannot have less points than eigenfunctions!');
    if sSimParams.b_randomStepSize
        vR = sort(randi([1 N],sParams.R,1));
    else
        step = N/sParams.R;
        vR = 1:step:N;
    end
    I = eye(sParams.M);
    
    mPhi_RM = mPhi(vR, :);
    vGR = vGi(vR);
    vCR = pinv(mPhi_RM) * vGR;
%     vCR = ( mPhi_RM.' * mPhi_RM + sParams.gamma*I ) \ (mPhi_RM.' * vGR);
%     vCR = vCR.*(abs(vCR) > 1e-15);
%     [(0:M-1)' vCR]
    vFi_hat = mPhi * vCR;
    accuracy = 100*(1 - norm(vFi_hat - vFi)/norm(vFi));
%     % ---------------------------------------------------------------------
%     % With kernel try
%     vFi_eval = vF_eval(:, i);
%     vPhi_m_x = zeros(M, 1);
%     vFi_hat2 = zeros(length(x_eval), 1);
%     for j = 1:length(x_eval)
%         for m = 0:M-1
%             vKernel_x_t = kernel(x_eval(j),t,l);
%             [vPhi_m_t, lambda_m] = phi(a, b, m, t);
%             vPhi_m_x(m+1) = trapz(t, vKernel_x_t.*vPhi_m_t.*vP_t);
%         end
%         vFi_hat2(j) = sum(vCR.*vPhi_m_x);
%         res = 10*exp(-2*x_eval(j)).*sin(5*x_eval(j));
%     end
%     % ---------------------------------------------------------------------
    cFigs{i} = figure(i+1);
    p1 = plot(x, vGi, 'Color',       '#4DBEEE', ...
                      'LineWidth',    2, ...
                      'DisplayName', ['$g(x)$']);
%     title({['$f_' num2str(i) '(x)$ with SNR = ' num2str(SNR) ' [dB]']; num2str(accuracy)}, 'Interpreter', 'latex', 'FontSize', 14)
    hold on
    p2 = plot(x, vFi, 'Color',       '#0072BD', ...
                      'LineWidth',    3, ...
                      'DisplayName', ['$f(x)$']);
    p3 = plot(x, vFi_hat, 'Color', '#D95319', ...
                          'LineWidth', 3, ...
                          'LineStyle', '-.', ...
                          'DisplayName', ['$\hat{f}(x)$']);

    p4 = plot(x(vR), vGi(vR), 'ro');
    hold off
    legend([p2 p1 p3], 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
    print(cFigs{i}, [sSimParams.outputFolder filesep 'fig' num2str(i+1) '_extrapolate_f' num2str(i)], '-depsc')
    fprintf('f%d : R = %d; M = %d; SNR = %.2f; Accuracy = %.2f%%\n', i, sParams.R, sParams.M, SNR, accuracy);
end
end