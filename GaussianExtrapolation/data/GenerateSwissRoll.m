function x = GenerateSwissRoll(N)

a = 1;   % swiss roll goes from a*pi to b*pi
b = 4;   
y = rand(2,N);
% uniform distribution along the manifold (in data space)
tt = sqrt((b*b-a*a)*y(1,:)+a*a);
tt = pi*tt;
% now tt should go from a*pi to b*pi
height = y(2,:);
x = [tt.*cos(tt)/b^2; height; tt.*sin(tt)/b^2]';


% r = round(sqrt(N));
% n = round(N/r); %round(sqrt(N));
% 
% t = (3*pi/2)*(1 + 2*linspace(0,1,n))';
% normFactor = [(2/(9*pi)), 1, (2/(9*pi))];
% x = normFactor.*[t.*cos(t), ones(n, 1), t.*sin(t)];
% xrep = repmat(x,r,1);
% 
% heights = [ones(n*r,1), kron(linspace(-1,1,r)',ones(n,1)), ones(n*r,1) ];
% x = xrep.*heights;
end