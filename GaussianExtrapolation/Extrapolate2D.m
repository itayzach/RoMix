function [] = Extrapolate2D(sParams, sSimParams)

assert(sParams.extrplM <= sParams.R, 'You cannot have less points than eigenfunctions!');
assert(sParams.dim == 2, 'This function works only for 2-D')

dx = 0.01;
x = (-2:dx:2-dx)';
N = length(x);

x1 = x.';
x2 = x.';
[mX1, mX2] = meshgrid(x1, x2);

%% F1
A1 = 3;
A4 = 0.001;
A6 = 0.0005;

phi1x1 = phi(sParams.a, sParams.b, 1, x1);
phi4x1 = phi(sParams.a, sParams.b, 4, x1);
phi6x1 = phi(sParams.a, sParams.b, 6, x1);

phi1x = phi1x1' * phi1x1;
phi4x = phi4x1' * phi4x1;
phi6x = phi6x1' * phi6x1;

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
    mPhi = zeros(N*N, sParams.extrplM);
    for m = 0:sParams.extrplM-1 
        vPhi_m_x1 = phi(sParams.a, sParams.b, m, x1);
        vPhi_m_x2 = phi(sParams.a, sParams.b, m, x2);
        
        % outter product since phi(x1,x2)=phi(x1)phi(x2)
        mPhi_m_x1x2 = vPhi_m_x1.' * vPhi_m_x2; 
        mPhi_m_x1x2 = mPhi_m_x1x2(:);
        mPhi(:, m+1) = mPhi_m_x1x2;
    end
    if sSimParams.b_randomStepSize
        vR = sort(randi([1 N^2],sParams.R,1));
    else
        step = floor(N^2/sParams.R);
        vR = 1:step:N^2;
                
    end
    I = eye(sParams.extrplM);
    
    mPhi_RM = mPhi(vR, :);
    mGi = reshape(vGi,N,N);
    vGR = vGi(vR);
    vCR = (mPhi_RM.' * mPhi_RM) \ ( mPhi_RM.' * vGR );
%     vCR = pinv(mPhi_RM) * vGR;

%     vC = zeros(sParams.extrplM, 1);
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
    fprintf('f%d : R = %d; M = %d; SNR = %.2f; Accuracy = %.2f%%\n', i, sParams.R, sParams.extrplM, SNR, accuracy);
end
end