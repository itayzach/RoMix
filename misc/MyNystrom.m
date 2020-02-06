%% Init
clc; clear; close all;
rng('default');
%% Read mnist data
% mnist = load('data/mnist.mat');
% X = single(mnist.testX);
% labels = mnist.testY.';
% X1 = X(labels == 1, :);
% X2 = X(labels == 2, :);
% data = [X1; X2];
% nPoints = size(data, 1);

%% Read corel data
% corel = load('data/corel_feature.mat');
% X = corel.feature;
% data = [X];
% nPoints = size(data, 1);

%% Create signal of time and space shifts
F = 1000;                           % Sampling frequency [Hz]
T = 1/F;                            % Sampling period    [sec]   
fCutoff = 10;                       % Cutoff frequency   [Hz]
d = 300;                            % Length of signal
t = (0:d-1);                        % Indexes vector
amp = 10;                           % Amplitude
signal = amp*sin(2*pi*fCutoff*t*T); % Signal
nPoints = 2000;                     % Number of signals

attenuation = exp(-linspace(0, 2, nPoints)).';%ones(nPoints, 1); %
attRep = repmat(attenuation, 1, length(signal));
sigRep = repmat(signal, length(attenuation), 1);
data = sigRep.*attRep + 0.1*amp*randn(size(sigRep)); 
figure; imagesc(data); colorbar; title('data'); 
%% Create a complete random data
% data = amp*randn(2000, 150);
%% Rand permutations
p = randperm(size(data, 1));
data = data(p, :);


%% Linear Kernel matrix
% W = data * data';
% W = ( W - min(min(W)) )/(max(max(W)) - min(min(W)));

%% Gaussian Kernel matrix
gaussSigma = 40;
% W = pdist2(data, data);
W = exp(-(pdist2(data, data).^2)./(2*gaussSigma^2));
% figure; imagesc((pdist2(data, data).^2)./(2*gaussSigma^2)); colorbar;
%% Verify W is Positive Definite
% figure; imagesc(W); colorbar; title('W');
e = eig(W);
fprintf('Max eigenvalue = %f\n', max(e)); 
fprintf('Min eigenvalue = %f\n', min(e));
fprintf('rank(W) = %d\n', rank(W));
assert(all(e > 0), 'W is not Positive Definite');

%% Set sizes
nSampledPointsVec = [500 ]; %300; %floor(1*r);
nOtherPointsVec = nPoints - nSampledPointsVec;
fprintf('nPoints = %d\n', nPoints);
fprintf('nSampledPointsVec = [ ');
fprintf('%d ', nSampledPointsVec);
fprintf(']\n');
%% Diffusion maps
D = diag(sum(W));

P = inv(D) * W;
L = D - W;

T = 5;
Pt = zeros(T, length(P), length(P));
Pt(1, :, :) = P;
for i = 2:T
    Pt(i, :, :) = squeeze(Pt(i-1, :, :)) * P;
end

% [Phi, Lambda, Psi] = svd(P);

% Psi are the right eigenvectors {P*Psi = Psi*Lambda   ==> P = Psi*Lambda*inv(Psi)}
% Lambda are the eigenvalues
% Phi are the left eigenvectors  {Phi'*P = Lambda*Phi' ==> P = inv(Phi')*Lambda*Phi'}

[Psi, Lambda, Phi] = eig(P);


%% Nystrom
normsVec = zeros(length(nSampledPointsVec), 1);
for i = 1:length(nSampledPointsVec)
    nSampledPoints = nSampledPointsVec(i);
    nOtherPoints = nOtherPointsVec(i);
    A = W(1:nSampledPoints, 1:nSampledPoints);
    B = W(1:nSampledPoints, nSampledPoints+1:end);
    C = B' * pinv(A) * B;
    Wh = [A B; B' C];
    normsVec(i) = norm(W - Wh);
%     figure; imagesc(abs(W - Wh)); colorbar; caxis([0 0.3]); title('|W(i,j) - Wh(i,j)|');
    
    d1 = sum([A; B'], 1);
    d2 = sum(B, 1) + sum(B, 2)'*pinv(A)*B;
    dhat = sqrt(1./[d1 d2])';
    A = A.*(dhat(1:nSampledPoints)*dhat(1:nSampledPoints)');
    B = B.*(dhat(1:nSampledPoints)*dhat(nSampledPoints+(1:nOtherPoints))');
    Asi = sqrtm(pinv(A));
    S = A + Asi*(B*B')*Asi;
    [U, L, T] = svd(S);
    V = [A; B']*Asi*U*pinv(sqrt(L));
    
%     E = zeros(size())
    nvec = 499;
    for j = 2:nvec+1
        E(:, j-1) = V(:, j)./V(:, 1);
    end
end
fprintf('normsVec = [ ');
fprintf('%f ', normsVec);
fprintf(']\n');

% Sort back ("undo" the random permutation)
Wsorted(p, :) = W;
WhSorted(p, :) = Wh;
figure; imagesc(Wsorted); colorbar; title('W');
figure; imagesc(WhSorted); colorbar; title('Wh');

keyboard