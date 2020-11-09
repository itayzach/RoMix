%% Implemetation of 
%  The nystrom extension for signals defined on a graph
%  Heimowitz, Ayelet, Eldar, Yonina C.
%% Restart
clc; clear; close all;
rng('default');

%% Parameters
b_showImages = false;

%% Load MNIST
sMNIST = load(fullfile('data', 'mnist.mat'));

% TODO: use the entire dataset
m = 5000;
trainX = sMNIST.trainX(1:m,:);
trainY = sMNIST.trainY(1:m).';
testX = sMNIST.testX(1:m,:);
testY = sMNIST.testY(1:m).';


[nTrain, dim] = size(trainX);
[nTest, ~] = size(testX);
N = nTest + nTrain;

xDim = sqrt(dim);
yDim = sqrt(dim);

%% Show some images
if b_showImages
    figure;
    colormap('gray')
    for i = 1:4
        randIdx = randi(nTrain);
        subplot(2,2,i);
        imagesc(reshape(trainX(randIdx,:), xDim, yDim).');
        title(['Train image #' num2str(randIdx) ' label = ' num2str(trainY(randIdx))]);
    end
end

%% Calculate distances
r = 1000;

k = 200;
F = squareform(pdist(trainX));
error('TODO: keep only the 200 lowest distances')
% [~, maxkIdx] = maxk(F, nTrain-k, 2);
% F(:, maxkIdx(2,:)) = 0;
W = nTrain^2 * F / sum(sum(F));
B = W(1:r, r+1:nTrain).';
d = sum(W,1);
De_sqrt = diag(sqrt(d(1:r)));
Db_msqrt = diag(1./sqrt(d(r+1:end)));


nDigits = 10;
% rIdx = randperm(nTrain, r);
% sR = zeros(r, nDigits);
sNystrom = zeros(nTrain-r, nDigits);
for i = 1:nDigits
    sR = (trainY(1:r) == i - 1);
    sNystrom(:,i) = Db_msqrt * B * De_sqrt * sR;
end

[~, argmax] = max(sNystrom, [], 2);
prediction = argmax - 1;
expected = trainY(r+1:end);

fprintf('Error: %.2f\n', sum(prediction ~= expected)/(nTrain-r))
