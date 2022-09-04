function [x, S] = GenerateSwissRoll(N, a, maxTheta, height, nGmmComp, b_plotSwissRoll, b_squeezeRollLayers)
% Running example:
%
% N = 5000; a = 1; maxTheta = 4*pi; height = 20; b_plotSwissRoll = true; nGmmComp = 5;
% [x, S] = GenerateSwissRoll(N, a, maxTheta, height, b_randn, b_plotSwissRoll);
%
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
if ~exist('nGmmComp', 'var')
    nGmmComp = 0;
end
if ~exist('b_squeezeRollLayers','var')
    b_squeezeRollLayers = false;
end

% construct archemedian spiral
% ============================
if b_squeezeRollLayers
    assert(maxTheta == 12*pi); % TD
end

% generate data
% ============================
% find angles which correspond to uniform sampling along spiral
if nGmmComp == 0
    theta = linspace(0, maxTheta, 100);
    s = SwissRollArclength(theta, a);
    sN = rand(N, 1)*max(s);
    thetaN = interp1(s, theta, sN);
else
    % split N equally among each component
    maxS = SwissRollArclength(maxTheta, a);
    %assert(nGmmComp == 5)
    sN = zeros(N,1);
    for c = 1:nGmmComp
        vInd = (1:N/nGmmComp) + (c-1)*(N/nGmmComp);
        sigma = maxS/(5*nGmmComp);
        mu = (c/(nGmmComp+1))*maxS;
        sN(vInd) = randn(N/nGmmComp, 1)*sigma + mu;
    end
    vRandInd = randperm(N);
    sN = sN(vRandInd);
    
    % transform
    theta = linspace(0, maxTheta, 100);
    s = SwissRollArclength(theta, a);
    thetaN = interp1(s, theta, sN);
    
    % verify
    assert(all(sN < maxS) && all(sN > 0))
end

% data uniformly distributed on swiss roll
if b_squeezeRollLayers
    x1 = a * cos(thetaN) .* thetaN.^0.5; %TD
    x2 = a * sin(thetaN) .* thetaN.^0.5; %TD
else
    x1 = a * cos(thetaN) .* thetaN;
    x2 = a * sin(thetaN) .* thetaN;
end
if nGmmComp > 0
    x3 = (height/5)*randn(N,1) + height/2;
else
    x3 = height*rand(N,1);
end

% store all data
if nGmmComp > 0
    S = [x1 x2 x3];
    x = [sN, x3];
else
    S = [sN, x3];
    x = [x1 x2 x3];
end
assert(all(~isnan(S(:))))
assert(all(~isnan(x(:))))

if b_plotSwissRoll
    figure; 
    subplot(131); plot(sN,thetaN,'o'); hold on; plot(sN,thetaN, '.'); xlabel('$s$'); ylabel('$\theta(s)$')
    subplot(132); scatter3(x1, x2, x3, 'filled'); xlabel('$x_1$'); ylabel('$x_2$'); zlabel('$x_3$'); view(30,80);
    subplot(133); scatter(sN, x3, 'filled'); xlabel('$\theta$'); ylabel('$ht$'); xlim([min(sN), max(sN)]); ylim([min(x3), max(x3)]);
    set(gcf,'Position', [250 400 1500 400])
end
end