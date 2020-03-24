function [] = Extrapolate2D(sParams, sSimParams)

assert(sParams.ExtrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
assert(sParams.dim == 2, 'This function works only for 2-D')

dx = 0.01;
x = (-2:dx:2-dx)';
N = length(x);

x1 = x.';
x2 = x.';
[mX1, mX2] = meshgrid(x1, x2);
X = [mX1(:) mX2(:)];
%% F1
A1 = 3;
A4 = 0.001;
A6 = 0.0005;

phi1x = reshape(phi(sParams, [0 1], X), length(x), length(x));
phi4x = reshape(phi(sParams, [2 0], X), length(x), length(x));
phi6x = reshape(phi(sParams, [3 1], X), length(x), length(x));

mF1      = A1*phi1x + A4*phi4x + A6*phi6x;
vF1      = mF1(:);
vF_awgn1 = sqrt(sSimParams.noiseVar1)*randn(N*N, 1);

%% F2
B1 = 5;
B2 = 30;
B3 = 30;

exp1 = B3*exp(-0.3*x1.^2).*sin(1*pi*x1);

mF2      = B2*exp(-0.5*mX1.^2).*exp(-0.5*mX2.^2).*sin(0.5*pi*mX1).*sin(0.5*pi*mX2); % + B3*exp(-0.3*mX1.^2).*exp(-0.3*mX2.^2).*sin(1*pi*mX1).*sin(1*pi*mX2);
vF2 = mF2(:);
vF_awgn2 = sqrt(sSimParams.noiseVar2)*randn(N*N, 1);

mF      = [vF1 vF2];
mF_awgn = [vF_awgn1 vF_awgn2];

%% Extrapolate
nFuncs  = size(mF, 2);
cFigs = cell(1, nFuncs);

for i = 1:nFuncs
    vFi = mF(:, i);
    vFi_awgn = mF_awgn(:, i);
    SNR = snr(vFi, vFi_awgn);
    vGi = vFi + vFi_awgn;
    mPhi = zeros(N*N, sParams.ExtrplM);
    for q = 0:sParams.ExtrplM-1 
        m = OneDim2TwoDimIndex(q, sParams.dim);
        vPhi_m_x = phi(sParams, m, X);
        mPhi(:, q+1) = vPhi_m_x;
    end
    if sSimParams.b_randomStepSize
        vR = sort(randi([1 N^2],sParams.R,1));
    else
        step = floor(N^2/sParams.R);
        vR = 1:step:N^2;
                
    end
    I = eye(sParams.ExtrplM);
    
    mPhi_RM = mPhi(vR, :);
    mGi = reshape(vGi,N,N);
    vGR = vGi(vR);
    vCR = (mPhi_RM.' * mPhi_RM) \ ( mPhi_RM.' * vGR );
%     vCR = pinv(mPhi_RM) * vGR;

%     vC = zeros(sParams.ExtrplM, 1);
%     vC(1+1) = A1;
%     fprintf('    vCR       vC \n')
%     disp([vCR, vC]);
    
        
    vFi_hat = mPhi * vCR;
    accuracy = 100*(1 - norm(vFi_hat - vFi)/norm(vFi));
    
    cFigs{i} = figure(i+1);
    subplot(1,2,1);
    p1 = surf(mX1, mX2, reshape(vGi,N,N), 'edgecolor', 'none');
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
    hold on;
    p4 = scatter3(mX1(vR), mX2(vR), mGi(vR), 'filled', 'ro');
    view(2)
%     view(90,0)
    title('$g({x})$', 'Interpreter', 'latex', 'FontSize', 14);
    colorbar;
    c1 = caxis;
    
    subplot(1,2,2);
    p2 = surf(mX1, mX2, reshape(vFi_hat,N,N), 'edgecolor', 'none');
    view(2)
%     view(90,0)
    title('$\hat{f}({x})$', 'Interpreter', 'latex', 'FontSize', 14);
    colorbar();
    xlabel('$x_1$', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$x_2$', 'Interpreter', 'latex', 'FontSize', 14);
    c2 = caxis;
    
    
    c3 = [min([c1 c2]), max([c1 c2])];
    caxis(c3)
    
    fig_left_loc = -1500;
    fig_bottom_loc = 400;
    fig_width = 1500;
    fig_height = 500;
    set(cFigs{i},'position',[fig_left_loc,fig_bottom_loc,fig_width,fig_height])
%     p1 = plot(x, vGi, 'Color',       '#4DBEEE', ...
%                       'LineWidth',    2, ...
%                       'DisplayName', ['$g(x)$']);
%     hold on
%     p2 = plot(x, vFi, 'Color',       '#0072BD', ...
%                       'LineWidth',    3, ...
%                       'DisplayName', ['$f(x)$']);
%     p3 = plot(x, vFi_hat, 'Color', '#D95319', ...
%                           'LineWidth', 3, ...
%                           'LineStyle', '-.', ...
%                           'DisplayName', ['$\hat{f}(x)$']);
% 
%     p4 = plot(x(vR), vGi(vR), 'ro');
%     hold off
%     legend([p2 p1 p3], 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
%     print(cFigs{i}, [sSimParams.outputFolder filesep 'fig' num2str(i+1) '_extrapolate_f' num2str(i)], '-depsc')
    fprintf('f%d : R = %d; M = %d; SNR = %.2f; Accuracy = %.2f%%\n', i, sParams.R, sParams.ExtrplM, SNR, accuracy);
end
end