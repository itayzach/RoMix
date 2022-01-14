function x = GenerateSwissRoll(N)

% According to:
% Parsimonious representation of nonlinear dynamical systems
% through manifold learning: A chemotaxis case study
% Carmeline J. Dsilva, 2015
% https://ronentalmon.com/wp-content/uploads/2019/03/ACHA_Dsilva_Jul_2015.pdf

% construct archemedian spiral
% ============================
a = 1;
% theta_vec = linspace(0, 4*pi, 100);
theta_vec = linspace(0, 12*pi, 100); %TD
s = 0.5*a*(theta_vec.*sqrt(1+theta_vec.^2)+log(theta_vec+sqrt(1+theta_vec.^2)));

% height
h = 20;

% generate data
% ============================
% find angles which correspond to uniform sampling along spiral
s_N = rand(N, 1)*max(s);
theta_N = interp1(s, theta_vec, s_N);

% data uniformly distributed on swiss roll
% x1 = a * cos(theta_N) .* theta_N;
% x2 = a * sin(theta_N) .* theta_N;
x1 = a * cos(theta_N) .* theta_N.^0.5; %TD
x2 = a * sin(theta_N) .* theta_N.^0.5; %TD
x3 = h*rand(N,1); 


% store all data
x = [x1 x2 x3]; 
% figure; plot(s,theta_vec,'o'); hold on; plot(s_N,theta_N, '.'); xlabel('s'); ylabel('\theta(s)')
% figure; scatter3(x1, x2, x3, 'filled'); xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); view(30,80);
end