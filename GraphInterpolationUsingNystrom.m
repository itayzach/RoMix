%% Init
clc; clear; close all;
rng('default');

%% Read mnist data
mnist = load('data/mnist.mat');
mXTest = single(mnist.testX);
% X_test = X_test(1:30,:);
mXTrain = single(mnist.trainX);
vTrainLabels = mnist.trainY.';
vTestLabels = mnist.testY.';
[nTestPoints, nPixels] = size(mXTest);
nTrainPoints = size(mXTrain, 1);
nPoints = nTestPoints + nTrainPoints;
N = nTestPoints;


%% Graph matricies
% F = pdist2([X_train; X_test], [X_train; X_test]);
mF = pdist2([mXTest], [mXTest]);
[mFSorted, mMN] = sort(mF, 2);
nSmallest = 200;
mMn = mMN(:,1:nSmallest);
mMnC = mMN(:, nSmallest+1:end);

sumF = 0;
for i = 1:N
   sumF = sumF + sum(mF(i, mMn(i, :)));
end

% -------------------------------------------------------------------------
% Calculate W from F
% -------------------------------------------------------------------------
mW = zeros(size(mF));
for i = 1:N
   mW(i,mMn(i, :)) = mF(i,mMn(i, :))*N^2 / sumF;
%    mW(i,mMn(i, :)) = 1-exp(-mF(i,mMn(i, :))*N^2 / sumF);
end

% -------------------------------------------------------------------------
% Make sure W is symmetric
% -------------------------------------------------------------------------
for i = 1:N
    for j = 1:N
        mW(i,j) = max(mW(i,j), mW(j,i));
    end
end

mD = diag(sum(mW));
%% Graph signals
nDigits = 10;
vS = zeros(N, nDigits);
for k = 0:9
    vS(:, k+1) = (vTestLabels == k);
end

%% Nystrom 
vR = [50 150 300 500 1000 1500 3000]; % number of sample points
vAccuracy = zeros(1, length(vR));
i = 1;
for r = vR
    mB = mW(r+1:N, 1:r);

    mDe = mD(1:r, 1:r);
    mDb = mD(r+1:N, r+1:N);


    vSSampled = vS(1:r, :);
    vSEst = diag(diag(mDb).^-0.5) * mB * diag(diag(mDe).^0.5) * vSSampled;
    [~, vSEstIdx] = max(vSEst, [], 2);

    vSampledTestLabels = vTestLabels(r+1:end);
    vSEstimatedLables = vSEstIdx - 1;

    vErrors = vSampledTestLabels ~= vSEstimatedLables;
    vAccuracy(i) = 100*(1-sum(vErrors)/N);

    fprintf('Accuracy = %.2f\n', vAccuracy(i))
    
    i = i + 1;
end

%% Plot accuracy graph
figure;
plot(vR, vAccuracy)
ylabel('Accuracy [%]')
xlabel('Size of sample set (r)')