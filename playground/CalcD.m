clc; clear; 

sSimParams = GetParameters();

dy = 0.01;
y = (-5:dy:5).';
nPoints = length(y);
D = zeros(nPoints, 1);
K_int = zeros(nPoints, 1);

assert(sSimParams.dim == 1)

for j = 1:nPoints
    integrand = @(x) kernel(sSimParams, y(j), x).*p(sSimParams, x);
    D(j) = integral(integrand,-1e3,1e3,'ArrayValued',true);
    
    integrand = @(x) kernel(sSimParams, y(j), x);
    K_int(j) = 1/sqrt(2*pi*sSimParams.ell^2) * integral(integrand,-1e3,1e3,'ArrayValued',true);
end
% D = (1/sqrt(2*pi*sSimParams.sigma))*D;

p_y = p(sSimParams, y);

figure;
plot(y, D, 'LineWidth', 2, 'DisplayName', '$\int K(x,y)p(y)dy$');
hold on
plot(y, p_y, 'LineWidth', 2, 'DisplayName', '$p(x)$');
plot(y, K_int, 'LineWidth', 2, 'DisplayName', '$\frac{1}{\sqrt{2 \pi \ell^2}} \int K(x,y)dy$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest')
title(['$\sigma =  ' num2str(sSimParams.sigma) '\quad \ell = ' num2str(sSimParams.ell) '$'], 'Interpreter', 'latex', 'FontSize', 14)
