function sTwoSpirals = GenerateTwoSpiralsDataset(nTrain, nTest)

noiseSigma = 0;
degrees = 570;
deg2rad = (2*pi)/360;
start = 90;
start = start * deg2rad;

[sTwoSpirals.x, sTwoSpirals.y] = GenerateTwoSpirals(nTrain, degrees, start, deg2rad, noiseSigma);
[sTwoSpirals.xt, sTwoSpirals.yt] = GenerateTwoSpirals(nTest, degrees, start, deg2rad, noiseSigma);

ell = 4;
posTrainInd = find(sTwoSpirals.y==1);
negTrainInd = find(sTwoSpirals.y==-1);
ipos = randperm(length(posTrainInd));
ineg = randperm(length(negTrainInd));
y1 = zeros(length(sTwoSpirals.y),1);
y1(posTrainInd(ipos(1:ell))) = 1;
y1(negTrainInd(ineg(1:ell))) = -1;
sTwoSpirals.y = y1;

function [data, labels] = GenerateTwoSpirals(N, degrees, start, deg2rad, noiseSigma)
% moon 1
N1 = floor(N/2);

n = start + sqrt(rand(N1,1)) * degrees * deg2rad;   
spiral_1x = -cos(n).*n + rand(N1,1)*noiseSigma;
spiral_1y = sin(n).*n + rand(N1,1)*noiseSigma;

% moon 2
N2 = N - N1;  
spiral_2x = cos(n).*n + rand(N2,1)*noiseSigma;
spiral_2y = -sin(n).*n + rand(N2,1)*noiseSigma;

% concate
data = 0.1*[spiral_1x ,spiral_1y; spiral_2x, spiral_2y];
labels = [ones(N1,1); -ones(N2,1)];

end

end

