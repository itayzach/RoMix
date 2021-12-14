% close all;
clear;
clc;

fid = fopen(fullfile('data', 'ml-100k', 'ua.base'));
fileCells = textscan(fid,'%d %d %d %d');
fclose(fid);

sData.userId    = fileCells{1};
sData.movieId   = fileCells{2};
sData.rating    = fileCells{3};
sData.timestamp = fileCells{4};


% sData = tdfread(fullfile('sData', 'ml-100k', 'ua.base'));
nUsers = max(sData.userId);
nMovies = max(sData.movieId);

R = NaN(nMovies,nUsers);
for row = 1:length(sData.userId)
    R(sData.movieId(row), sData.userId(row)) = sData.rating(row);
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

% sparse_matrix = R.';
% min_common = 2;
% adjacency = zeros(nMovies, nMovies);
% for i = 1:nMovies
%     for j = 1:i-1
%         sparse_i = R(i,:);
%         sparse_j = R(j,:);
%         common_idx = ~isnan(sparse_i .* sparse_j);
%         if sum(common_idx) < min_common
%             adjacency(i,j) = 0;
%             adjacency(j,i) = 0;
%             continue
%         end
%         sparse_i(~common_idx) = NaN;
%         sparse_j(~common_idx) = NaN;
%         sparse_i(isnan(sparse_i)) = 0;
%         sparse_j(isnan(sparse_j)) = 0;
%         adjacency(i, j) = sum(sparse_i .* sparse_j / norm(sparse_i) / norm(sparse_j));
%         adjacency(j, i) = adjacency(i, j);
%     end
%     if (mod(i,100) == 0)
%         fprintf('iteration %d/%d\n', i, nMovies)
%     end
% end
% % adjacency(isnan(adjacency)) = 0;
% degree_matrix = diag(sum(adjacency, 1));
% laplacian = degree_matrix - adjacency;
% 
% predict(sparse_matrix, laplacian)
% 
% 
% % function sPrediction = predict(sparse_matrix, laplacian)
% idx = find(diag(laplacian) ~= 0);
% laplacian = laplacian(idx,:);
% laplacian = laplacian(:,idx);
% predicted_matrix = zeros(size(sparse_matrix));
% theta = nanmean(sparse_matrix);
% degree_matrix = diag(diag(laplacian));
% degree_matrix_half = diag(diag(degree_matrix) .^ (-1/2));
% degree_matrix_half(degree_matrix_half==inf) = 0;
% n_laplacian = degree_matrix_half * laplacian * degree_matrix_half;
% assert(~any(isnan(n_laplacian(:))))
% laplacian_square = n_laplacian*n_laplacian;
% 
% figure; plot(n_laplacian(:,1:10),'.')

% for i = 1:
% du_signal = sparse_matrix(1,idx);
% known_idx = ~np.isnan(du_signal)
% unknown_idx = np.isnan(du_signal)
% du_signal[np.isnan(du_signal)] = 0
% laplacian_sc = laplacian_square[~known_idx, :]
% laplacian_sc = laplacian_sc[:, ~known_idx]
% e_value, _ = np.linalg.eig(laplacian_sc)
% sorted_indices = np.argsort(e_value)
% W = np.sqrt(e_value[sorted_indices[0]])
% % Step 3
% 
% 
% 
% predicted_matrix(predicted_matrix<0) = 0;
% predicted_matrix(predicted_matrix>5) = 5;
% sPrediction.predicted_matrix = predicted_matrix;
% sPrediction.idx = idx;
% end



