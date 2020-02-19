a = 4;
b = 4.08;

x1 = -0.5:0.05:0.5;
x2 = x1;
[XX1,XX2] = meshgrid(x1,x2);
M = 4;

figure;
for m = 0:M-1 
    [vPhi_m_x1, lambda_m1] = SqExpEig(a, b, m, x1);
    [vPhi_m_x2, lambda_m2] = SqExpEig(a, b, m, x2);
    
    F = vPhi_m_x1.' * vPhi_m_x2; 
    subplot(2,2,m+1);
    surf(XX1,XX2,F)
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$\phi(x_1,x_2)$', 'Interpreter', 'latex')
end