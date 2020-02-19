function [f,labels,error]=ml_test_eigrls(classifier,XTest, YTest)

% ML_TEST Uses a classifier to classify data in XTest
% ----------------------------------------------------------------------------%
%
% Usage: 
% [f,labels,error]=ml_test(classifier,XTest,YTest)
%
% Inputs:
% classifier: A classifier structure returned by ml_train or saveclassifier
% XTest : n XTest d matrix (n examples in d dimensions)
% YTest : optional column vector of labels. Values in -1,0,+1 
%     If YTest is provided, computes error rates 
%     Note:  YTest is allowd to be in [-1,0,+1]
%            The error computation is done over labeled points [-1,+1]
% 
% Outputs: 
% f : real valued classifier output
% labels: f thresholded at b 
% error: error rate over labeled part of YTest
% 
% Author:  Vikas Sindhwani vikass@cs.uchicago.edu
%          June 2004
%------------------------------------------------------------------------------%


% read classifier
a_k=classifier.a_k;
b_k=classifier.b_k;
alpha=classifier.alpha;
M=classifier.M;
xtrain=classifier.xtrain;
ytrain=classifier.ytrain;

for m = 0:M-1 
    [vPhi_m_x1, lambda_m1] = SqExpEig(a_k, b_k, m, xtrain(:,1));
    [vPhi_m_x2, lambda_m2] = SqExpEig(a_k, b_k, m, xtrain(:,2));
    
    Phi(:,m+1) = vPhi_m_x1 .* vPhi_m_x2; 
%     subplot(2,2,m+1);
%     surf(XX1,XX2,F)
%     xlabel('$x_1$', 'Interpreter', 'latex')
%     ylabel('$x_2$', 'Interpreter', 'latex')
%     zlabel('$\phi(x_1,x_2)$', 'Interpreter', 'latex')
end
f=Phi*alpha;
labels=sign(f);

% compute error rate over labeled part of test set
if exist('YTest','var')==1
    test=find(YTest);
    error=sum(labels(test)~=YTest(test))/length(test)*100;
end

