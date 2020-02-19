function plotclassifiers(classifier,xrange, yrange);

x=xrange;
y=yrange;

%x=-1.5:0.05:2.5; % x range
%y=-1.5:0.05:2.5; % y range

%x=0:0.01:1; % x range
%y=0:0.01:1; % y range


KERNELTYPE={'linear','poly','rbf'};

alpha=classifier.alpha;
Xtrain=classifier.xtrain;
b=classifier.b;
method=classifier.Name;
kerneltype=classifier.Kernel;
kernelparam=classifier.KernelParam;
lambda=classifier.gammas;


[xx,yy]=meshgrid(x,y); % makes a grid of points

X=[xx(:) yy(:)];


if (strcmp(method,'linear_rlsc') | strcmp(method,'linear_laprlsc'))

    z=X*alpha + b;

elseif strcmp(method,'eigrls')
    
    for m = 0:classifier.M-1 
        [vPhi_m_x1, ~] = SqExpEig(classifier.a_k, classifier.b_k, m, X(:,1));

        [vPhi_m_x2, ~] = SqExpEig(classifier.a_k, classifier.b_k, m, X(:,2));

%         lambda(m+1) = lambda_m1;
        Phi(:, m+1) = vPhi_m_x1.*vPhi_m_x2;           
    end
    
    z=Phi*alpha;
    
    z = sign(z);
else
    K=calckernel(kerneltype,kernelparam, Xtrain, X);
    
    Ka = K*alpha;
    z=sign(Ka-b);
    
 
end

Z=reshape(z,length(x),length(y));


[cs,h]=contourf(xx,yy,Z,[0 0]);shading flat;
%clabel(cs,h);colorbar

hold on;
title([method '  ('  kerneltype  ') : ' '\sigma =' num2str(kernelparam) ], 'FontSize', 14);
if size(lambda,2)==2
  xlabel([' \gamma_A = ' sprintf('%.3f', lambda(1)) '      \gamma_I = ' num2str(lambda(2))], 'FontSize', 14);
else
  xlabel([' \gamma_I = ' sprintf('%.3f', lambda(1))], 'FontSize', 14);  
end

% surf(xx,yy,Z)
% hold on;
% xlabel('$x_1$', 'Interpreter', 'latex')
% ylabel('$x_2$', 'Interpreter', 'latex')
% zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
% colorbar;
%xlabel('White (-1) Green (+1)');
%title( [' \lambda_1 = ' num2str(lambda(1)) '   \lambda_2 = '  num2str(lambda(2))]);

%xlabel([' \lambda_1 = ' num2str(lambda(1)) ]);
%ylabel(['\lambda_2 = '  num2str(lambda(2))]);
 
%plot2D(Xtest,Ytest,5);
