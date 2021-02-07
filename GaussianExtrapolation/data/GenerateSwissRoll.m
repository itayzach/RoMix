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

end