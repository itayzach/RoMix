clear; close all; clc;
%% Random data
% n = 2000; % Number of time points
% d = 20000; % Dimension
% k = 10000; % Random dimension
% 
% X = 10+5*randn(n,d); % Data matrix


%% MNIST
mnist = load('data/mnist.mat');
X = single(mnist.testX);
Y = mnist.testY.';
zeroInd = find(Y == 1);
oneInd = find(Y == 5);
% twoInd = find(Y == 2);
ind = [zeroInd; oneInd];
X = X(ind,:);
Y = Y(ind);
[n, d] = size(X);

figure;
colormap('gray')
for i = 1:4
    subplot(2,2,i);
    imagesc(reshape(X(980+i,:), sqrt(d), sqrt(d)).');
    title(['Train image #' num2str(980+i) ' label = ' num2str(Y(980+i))]);
end
%% RP
k = 3;
R = randn(d,k); % Random projection matrix
for kkk = 1:d
   R(kkk,:) = R(kkk,:)/norm(R(kkk,:)); % Normalization
end

P = X*R; % Projection

% Distances
ii = 10; jj = 20;
norm(X(ii,:)-X(jj,:))

norm(P(ii,:)-P(jj,:))

figure;
c = Y;
scatter3(P(:,1),P(:,2),P(:,3),[],c,'filled')

figure;
subplot(2,1,1); histogram(X(:),100); title('orig histogram', 'Interpreter', 'latex', 'FontSize', 14); set(gca,'FontSize', 14);
subplot(2,1,2); histfit(P(:),100); title('new histogram', 'Interpreter', 'latex', 'FontSize', 14); set(gca,'FontSize', 14);
set(gcf,'Position', [800 100 400 600])