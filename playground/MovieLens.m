close all;
clear;
clc;


data = tdfread(fullfile('data', 'ml-100k', 'u1.base'));
nUsers = max(data.userId);
nMovies = max(data.movieId);

R = zeros(nMovies,nUsers);
for row = 1:length(data.userId)
    R(data.movieId(row), data.userId(row)) = data.rating(row);
end

% remove movies that have no ratings at all
moviesToRemove = sum(R,2)==0;
R(moviesToRemove,:) = [];
nMovies = size(R,1);


% cosine similarity
S = R*R.';
figure; imagesc(S); colorbar;
Q = diag(1./vecnorm(R,2,2));
W = Q*S*Q;

% exp
% pdistR = pdist(R);
% omega = sqrt(median(pdistR));
% dist = squareform(pdistR);
% W = exp(-dist/(2*omega^2));

% Y = tsne(R,'Algorithm','exact','Distance','cosine');
% figure; gscatter(Y(:,1),Y(:,2))


% laplacian
d = sum(W,2);
% L = eye(nMovies) - diag(d.^(-1/2))*W*diag(d.^(-1/2));
L = diag(d) - W;


assert(~any(isnan(S(:))))
[Psi,Lambda] = eig(L);
lambda = diag(Lambda);
figure; plot(real(Psi(:,1:5)));
