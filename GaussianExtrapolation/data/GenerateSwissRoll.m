function x = GenerateSwissRoll(N)

%% Method #1
% a = 1;   % swiss roll goes from a*pi to b*pi
% b = 4;   
% y = rand(2,N);
% % uniform distribution along the manifold (in data space)
% tt = sqrt((b*b-a*a)*y(1,:)+a*a);
% tt = pi*tt;
% % now tt should go from a*pi to b*pi
% height = y(2,:);
% x = [tt.*cos(tt)/b^2; height; tt.*sin(tt)/b^2]';


%% Method #2
% r = round(sqrt(N));
% n = round(N/r); %round(sqrt(N));
% 
% t = (3*pi/2)*(1 + 2*linspace(0,1,n))';
% normFactor = [(2/(9*pi)), 1, (2/(9*pi))];
% x = normFactor.*[t.*cos(t), ones(n, 1), t.*sin(t)];
% xrep = repmat(x,r,1);
% 
% heights = [ones(n*r,1), kron(linspace(-1,1,r)',ones(n,1)), ones(n*r,1) ];
% x = xrep.*heights;



%% Method #3
% According to:
% Parsimonious representation of nonlinear dynamical systems
% through manifold learning: A chemotaxis case study
% Carmeline J. Dsilva, 2015
% https://ronentalmon.com/wp-content/uploads/2019/03/ACHA_Dsilva_Jul_2015.pdf
% construct archemedian spiral
% ============================
a = 1;
theta_vec = linspace(0, 4*pi, 100);
s = 0.5*a*(theta_vec.*sqrt(1+theta_vec.^2)+log(theta_vec+sqrt(1+theta_vec.^2)));

% height
h = 20;

% generate data
% ============================
% % intialize random number generator
% rng(321);

% find angles which correspond to uniform sampling along spiral
theta = interp1(s, theta_vec, rand(N, 1)*max(s));

% data uniformly distributed on swiss roll
x1 = a * cos(theta) .* theta;
x2 = a * sin(theta) .* theta;
x3 = h*rand(N,1); 


% store all data
x = [x1 x2 x3]; 
end