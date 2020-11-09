clc; close all; clear;

z = -10:0.1:10;
J = zeros(5,201);
for i = 0:4
    J(i+1,:) = real(sqrt(pi./(2*z)).*besselj(i+0.5,z));
    J(i+1,z < 0) = -J(i+1,z < 0);
    if i == 0
        J(i+1, z == 0) = 1;
    else
        J(i+1, z == 0) = 0;
    end
end
plot(z,J)
grid on
legend('J_0','J_1','J_2','J_3','J_4','Location','Best')
title('Bessel Functions of the First Kind for $\nu \in [0, 4]$','interpreter','latex')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')