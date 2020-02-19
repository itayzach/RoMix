function [alpha, XTrain]=experiment_moon(XTrain,YTrain,XTest,YTest,method,q,s, options);

% 2 Moons Experiment
% Author: Vikas Sindhwani (vikass@cs.uchicago.edu)




%q=-1:5;
if nargin==6 % perform a search over lambda1 and lambda2
    lambda1=2.^q;
    lambda2=2.^q;
else % interpret q as lambda1 and s as lambda2
    
    lambda1=q;
    lambda2=s;
    
end

% best SVM
min_err=100;

if strcmp(method,'svm') || strcmp(method,'tsvm') || strcmp(method,'rlsc')
    % optimize over just 1 parameter
    for i=1:length(q)
        options.gamma_A=lambda1(i);
        options.gamma_I=0;
        classifier=ml_train(XTrain,YTrain,options, method);
        [f,labels,error]=ml_test(classifier,XTest,YTest);
        fprintf('lambda(i) = %.4f, error = %.4f, min_err = %.4f\n', lambda1(i), error, min_err);
        if error < min_err
            min_err=error;
            best_classifier=classifier;
        end
        
    end
    
    
elseif strcmp(method, 'eigrls')
    for i=1:length(q)
        for j=1:length(q)
            options.gamma_A=lambda1(i);
            options.gamma_I=lambda2(j);
            classifier=ml_train(XTrain,YTrain,options, method);
            [f,labels,error]=ml_test_eigrls(classifier,XTest,YTest);
            fprintf('lambda(i) = %.4f, error = %.4f, min_err = %.4f\n', lambda1(i), error, min_err);
            if error <= min_err
                min_err=error;
                best_classifier=classifier;
            end
        end
    end
else % optimize over both
    
    for i=1:length(q)
        for j=1:length(q)
            options.gamma_A=lambda1(i);
            options.gamma_I=lambda2(j);
            classifier=ml_train(XTrain,YTrain,options, method);
            [f,labels,error]=ml_test(classifier,XTest,YTest);
            fprintf('lambda(i) = %.4f, error = %.4f, min_err = %.4f\n', lambda1(i), error, min_err);
            if error <= min_err
                min_err=error;
                best_classifier=classifier;
            end
        end
    end
    
end
lab=find(YTrain);
xmin=min(XTrain(:,1)); ymin=min(XTrain(:,2)); rmin=min(xmin,ymin)-0.2;
xmax=max(XTrain(:,1)); ymax=max(XTrain(:,2)); rmax=max(xmax,ymax)+0.2;
steps=(rmax-rmin)/100;
xrange=rmin:steps:rmax;
yrange=rmin:steps:rmax;
plotclassifiers(best_classifier, xrange, yrange);

alpha = best_classifier.alpha;
XTrain = best_classifier.xtrain;






% i=0; j=0;
%
% for lgl2=1:nq
%  i=i+1; j=0;
%
%     for lgl1=1:nq
%  j=j+1;
%  p=p+1;
% options.lambda1=10^(q(lgl1));
% options.lambda2=10^(q(lgl2));
% classifier=Train(XTrain,Y1,options,method);
% [fT,labels,error_rate]=Test(classifier, XTest,YTest);
% subplot(nq,nq,p); plotclassifiers(XTest,YTest,classifier);plot2D(XTrain(lab,:),YTrain(lab),6);
% e(i,j)=error_rate;
%
% end
% end
