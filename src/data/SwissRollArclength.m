function s = SwissRollArclength(theta, a)
if ~exist('a', 'var')
    a = 1;
end
s = 0.5*a*(theta.*sqrt(1+theta.^2)+log(theta+sqrt(1+theta.^2)));
end