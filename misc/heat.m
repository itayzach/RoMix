%Delete all stored variables and clear the command window. 
clear
clc
%PARAMETERS:
%Number of terms in Fourier series
N=500;
%Thermal diffusivity
diffusivity=0.000111;
% x axis.
x=linspace(0,1,100);
% Placeholder for graph data
y=zeros(1,length(x));
%These will be used to fix the graph
X=linspace(-0.15,1.15,10);
Y=linspace(0,120,12);
%Time axis. 
t=linspace(0,60,60);
%Solution matrix placeholder. Each row will represent T(x,t) for a value of t.
soln=zeros(length(t), length(y));
%Fourier coefficients placeholder.
coeffs = zeros(1,N);
%Temperature data
data=[0 0 0 0 0 1 0 0 0 0 ];

%COMPUTATION:

%Compute the coefficients.
for n=1:N
        for k=1:10
            zmin=0.1*k-0.1;
            zmax=0.1*k;
            A=data(k);
            coeffs(n) = coeffs(n)+integral(@(z)2*A*sin(n*pi*z),zmin,zmax);
        end
        y = y+coeffs(n)*sin(n*pi*x);
end

% Compute the solution.
for a=1:length(t)
    T=t(a);
    y=zeros(1,length(x));
    for n=1:N
        y = y+coeffs(n)*sin(n*pi*x)*exp(-1*diffusivity*(n^2)*(pi^2)*T);
    end
    %The temperature should be greater than or equal to zero and less than or equal to 100, this fixes
    %any parts of the data that are slightly out of range due to approximation
    %error.
    for k = 1:length(y)
        if (y(k)<0);
            y(k) = 0;
        elseif (y(k)>100);
            y(k)=100;
        else
            y(k)=y(k);
        end
    end
    %Row number a of solution matrix is equal to the y array.
    soln(a,:)=flip(y);
end
figure;
mesh(x,t,soln)
xlabel('x');
ylabel('t');
zlabel('T(x,t) in Celsius');
xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
xticklabels({'1', '0.8', '0.6', '0.4', '0.2', '0'});
%This part produces the .gif animation. I don't know for sure if this part will work in
%Scilab without modification.
% fig=figure;
% elapsed=0;
% nImages=length(t);
% for T = 1:nImages
%     clf
%     if (t(T)-(elapsed+1))>=0
%         elapsed = elapsed+1;
%     end
%     FF=zeros(1,length(Y))-0.125;
%     hold on
%     plot(x,soln(T,:))
%     plot(X,zeros(1,length(X)), 'black');
%     plot(FF,Y, 'white');
%     hold off
%     title(['Elapsed time = ' num2str(elapsed) ' seconds'])
%     frame = getframe(fig);
%     im{T} = frame2im(frame);
% end
% close
% filename = 'RandomHeatSim1.gif';
% for idx = 1:nImages
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1;
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1.5);
%     elseif (idx > 1) && (idx < nImages);
%         if (idx>200) && (mod(idx,3)==0) && (idx <= 401) %After 200 frames, 3x frameskip.
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
%         elseif (idx<=200)
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
%         elseif (idx>401) && (mod(idx,5)==0)  && (idx <900) %After 401 frames, 5x frameskip.
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
%         elseif (idx>=901) && (mod(idx,10)==0) %After 900 frames, 10x frameskip.
%             imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.02);
%         end
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',1.5);
%     end
% end