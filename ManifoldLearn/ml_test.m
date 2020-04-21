function [f,labels,error]=ml_test(classifier,X,Y)

% ML_TEST Uses a classifier to classify data in X
% ----------------------------------------------------------------------------%
%
% Usage: 
% [f,labels,error]=ml_test(classifier,X,Y)
%
% Inputs:
% classifier: A classifier structure returned by ml_train or saveclassifier
% X : n x d matrix (n examples in d dimensions)
% Y : optional column vector of labels. Values in -1,0,+1 
%     If Y is provided, computes error rates 
%     Note:  Y is allowd to be in [-1,0,+1]
%            The error computation is done over labeled points [-1,+1]
% 
% Outputs: 
% f : real valued classifier output
% labels: f thresholded at b 
% error: error rate over labeled part of Y
% 
% Author:  Vikas Sindhwani vikass@cs.uchicago.edu
%          June 2004
%------------------------------------------------------------------------------%

alpha = classifier.alpha;
xtrain = classifier.xtrain;

if strcmp(classifier.Name, 'eigrls')
    sParams = GetParameters();

    Phi_xtrain = zeros(size(xtrain,1), sParams.ExtrplM);
    Phi_xtest = zeros(size(X,1), sParams.ExtrplM);
    lambda_m = zeros(sParams.ExtrplM, 1);
    for i = 0:sParams.ExtrplM-1 
        m = OneDim2TwoDimIndex(i);
        lambda_m(i+1) = lambda(sParams, m);
        Phi_xtrain(:,i+1) = phi(sParams, m, xtrain);
        Phi_xtest(:,i+1) = phi(sParams, m, X);
    end
    PLP = Phi_xtest*diag(lambda_m)*Phi_xtrain.';
    f_PLP_alpha=PLP*alpha;

    c = classifier.c;
    c_from_alpha = diag(lambda_m)*Phi_xtrain.'*alpha;
    isalmostequal(c,c_from_alpha,1e-10,'',false);
    f_Pc_from_alpha = Phi_xtest*c_from_alpha;
    f_Pc = Phi_xtest*c;
    
    Kernel=classifier.Kernel;
    KernelParam=classifier.KernelParam;
    K=calckernel(Kernel,KernelParam,xtrain,X);
    f_alpha = K*alpha;
    
    isalmostequal(K,PLP,1e-10,'',false);
    
    
%     figure;
%     subplot(2,1,1);
%         imagesc(K); colorbar;
%         title('kernel');
%     
%     subplot(2,1,2);
%         imagesc(PLP); colorbar;
%         title('mercer');
    f = f_Pc;
    isalmostequal(f,f_alpha,1e-10,'',false);
    
else
    % read classifier
    Kernel=classifier.Kernel;
    KernelParam=classifier.KernelParam;


    % compute test kernel
    K=calckernel(Kernel,KernelParam,xtrain,X);
    f=K*alpha;
end
labels=sign(f);

% compute error rate over labeled part of test set
if exist('Y','var')==1
    test=find(Y);
    error=sum(labels(test)~=Y(test))/length(test)*100;
end

