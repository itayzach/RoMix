%% load data

load('mnist.mat');
[l,k] = size(testX);
m = 10;
% [k,l,m] = size(data); % k = 256, l = 1100, m = 10

for iter = 1:10
    
    fprintf('Generating data-set %d...\n', iter);
    
    % number of samples to consider from each class
    num_samples = 100;
    N = m*num_samples;

    % graph signal
    f = zeros(N,1);
    X = zeros(N,k);
    for i = 1:m
        classInd = find(testY == i-1, num_samples);
        f((i-1)*num_samples + 1 : (i-1)*num_samples + num_samples) = i;
        X((i-1)*num_samples + 1 : (i-1)*num_samples + num_samples, :) = testX(classInd, :);
%         figure; imagesc(reshape(X((i-1)*num_samples + 1,:),28,28).')
    end

    % membership functions
    mem_fn = false(N,m);
    for i = 1:m
        mem_fn(f==i,i) = true;
    end
    
    % save
    save(['mnist_set' num2str(iter) '.mat'],'mem_fn', 'X');

end
