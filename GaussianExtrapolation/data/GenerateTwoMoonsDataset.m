function sTwoMoons = GenerateTwoMoonsDataset(nTrain, nTest, noiseSigma)

if ~exist('sigmad', 'var')
    noiseSigma = 0.1;
end

[sTwoMoons.x, sTwoMoons.y] = GenerateTwoMoons(nTrain, noiseSigma);
[sTwoMoons.xt, sTwoMoons.yt] = GenerateTwoMoons(nTest, noiseSigma);

ell = 1;
pos = find(sTwoMoons.y==1);
neg = find(sTwoMoons.y==-1);
ipos = randperm(length(pos));
ineg = randperm(length(neg));
y1 = zeros(length(sTwoMoons.y),1);
y1(pos(ipos(1:ell)))=1;
y1(neg(ineg(1:ell)))=-1;
sTwoMoons.y = y1;


function [data, labels] = GenerateTwoMoons(n, noiseSigma)
% moon 1
N1 = floor(n/2);
moon1x = cos(linspace(0, pi, N1)) + noiseSigma*randn(1,N1);
moon1y = sin(linspace(0, pi, N1)) + noiseSigma*randn(1,N1);

% moon 2
N2 = n - N1;
moon2x = 1 - cos(linspace(0, pi, N2)) + noiseSigma*randn(1,N2);
moon2y = 1 - sin(linspace(0, pi, N2)) + noiseSigma*randn(1,N2) - 1;

% concate
data = [moon1x.',moon1y.'; moon2x.', moon2y.'];
labels = [-ones(N1,1); ones(N2,1)];

end

end

