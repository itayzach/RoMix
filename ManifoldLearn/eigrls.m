% function [alpha,b]=eigrls(Y,Phi, Lambda,gamma1,gamma2)
function [alpha, b]=eigrls(L, Y, Phi, Lambda, gamma1,gamma2)

n=size(Y,1); % total examples
l=sum(abs(Y)); % labeled examples
M = size(Phi,2);

u = n - l; % unlabeled examples
I=eye(M);
J=diag(abs(Y));


if ~isempty(Lambda) && gamma2~=0
    alpha = ( (J*Phi).' * (J*Phi) + gamma1*I + gamma2*Lambda ) \ ((J*Phi).' * Y);
else
    alpha = ( (J*Phi).' * J*Phi + gamma1*I ) \ ((J*Phi).' * Y);
end
b = 0;
end
