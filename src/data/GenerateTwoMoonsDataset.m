function sTwoMoons = GenerateTwoMoonsDataset(nTrain, nTest, nLabeled, b_loadTwoMoonsMatFile, noiseSigma)

if ~exist('noiseSigma', 'var')
    noiseSigma = 0.1;
end

if b_loadTwoMoonsMatFile
    sTwoMoons = load(strcat('data', filesep, '2moons.mat'));
else
    [sTwoMoons.x, sTwoMoons.y] = GenerateTwoMoons(nTrain, noiseSigma);
    [sTwoMoons.xt, sTwoMoons.yt] = GenerateTwoMoons(nTest, noiseSigma);
end

ell = nLabeled/2;
posTrainInd = find(sTwoMoons.y==1);
negTrainInd = find(sTwoMoons.y==-1);
ipos = randperm(length(posTrainInd));
ineg = randperm(length(negTrainInd));
y1 = zeros(length(sTwoMoons.y),1);
y1(posTrainInd(ipos(1:ell))) = 1;
y1(negTrainInd(ineg(1:ell))) = -1;
sTwoMoons.y = y1;

b_driftMoonsAway = true;
if b_loadTwoMoonsMatFile && b_driftMoonsAway
    % Increase the distance between the two moons
    posTestInd = find(sTwoMoons.yt==1);
    negTestInd = find(sTwoMoons.yt==-1);
    sTwoMoons.xt(posTestInd,2) = sTwoMoons.xt(posTestInd,2) + 0.25;
    sTwoMoons.xt(negTestInd,2) = sTwoMoons.xt(negTestInd,2) - 0.25;
    sTwoMoons.x(posTrainInd,2) = sTwoMoons.x(posTrainInd,2) + 0.25;
    sTwoMoons.x(negTrainInd,2) = sTwoMoons.x(negTrainInd,2) - 0.25;
end


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
labels = [ones(N1,1); -ones(N2,1)];

end

end

