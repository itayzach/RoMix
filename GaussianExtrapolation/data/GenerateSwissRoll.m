function [x, S] = GenerateSwissRoll(N, a, maxTheta, height, b_plotSwissRoll, b_squeezeRollLayers)
% According to:
% Parsimonious representation of nonlinear dynamical systems
% through manifold learning: A chemotaxis case study
% Carmeline J. Dsilva, 2015
% https://ronentalmon.com/wp-content/uploads/2019/03/ACHA_Dsilva_Jul_2015.pdf

if ~exist('a', 'var')
    a = 1;
end
if ~exist('maxTheta', 'var')
    maxTheta = 4*pi;
end
if ~exist('height', 'var')
    height = 20;
end
if ~exist('b_plotSwissRoll', 'var')
    b_plotSwissRoll = false;
end

if ~exist('b_squeezeRollLayers','var')
    b_squeezeRollLayers = false;
end

% construct archemedian spiral
% ============================
if b_squeezeRollLayers
    assert(maxTheta == 12*pi); % TD
end
theta = linspace(0, maxTheta, 100);
s = SwissRollArclength(theta, a);

% generate data
% ============================
% find angles which correspond to uniform sampling along spiral
sN = rand(N, 1)*max(s);
thetaN = interp1(s, theta, sN);

% data uniformly distributed on swiss roll
if b_squeezeRollLayers
    x1 = a * cos(thetaN) .* thetaN.^0.5; %TD
    x2 = a * sin(thetaN) .* thetaN.^0.5; %TD
else
    x1 = a * cos(thetaN) .* thetaN;
    x2 = a * sin(thetaN) .* thetaN;
end
x3 = height*rand(N,1);

% store all data
S = [sN, x3];
x = [x1 x2 x3];

if b_plotSwissRoll
    figure; 
    subplot(121); plot(sN,thetaN,'o'); hold on; plot(sN,thetaN, '.'); xlabel('s'); ylabel('\theta(s)')
    subplot(122); scatter3(x1, x2, x3, 'filled'); xlabel('x_1'); ylabel('x_2'); zlabel('x_3'); view(30,80);
    set(gcf,'Position', [450 400 1000 400])
end
end