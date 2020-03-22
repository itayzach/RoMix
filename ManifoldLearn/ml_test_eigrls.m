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
% xtrain=classifier.xtrain;
% ytrain=classifier.ytrain;

sParams.dim = 2;
sParams.constsType = 1;
sParams.a = a_k;
sParams.b = b_k;
sParams.ell = 1/sqrt(2*sParams.b); % kernel width
sParams.sigma = 1./(2*sParams.a);
sParams.mu = 0*ones(1, sParams.dim);
sParams.c = sqrt(sParams.a.^2 + 2*sParams.a.*sParams.b);
sParams.A = sParams.a + sParams.b + sParams.c;
sParams.B = sParams.b./sParams.A;

Phi = zeros(size(XTest,1), M);
for i = 0:M-1 
    m = OneDim2TwoDimIndex(i, sParams.dim);
    Phi(:,i+1) = phi(sParams, m, XTest);
end
f=Phi*alpha;
labels=sign(f);

% compute error rate over labeled part of test set
if exist('YTest','var')==1
    test=find(YTest);
    error=sum(labels(test)~=YTest(test))/length(test)*100;
end

