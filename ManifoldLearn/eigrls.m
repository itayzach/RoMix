function [alpha, b]=eigrls(Y, Phi, Lambda, gamma1, gamma2, L)

n = size(Y,1);   % total examples
l = sum(abs(Y)); % labeled examples
u = n - l;       % unlabeled examples

M = size(Phi,2);
% I = eye(M);
I = eye(n);
J = diag(abs(Y));

% % alpha = ( (J*Phi).' * (J*Phi) + gamma1*I + gamma2*Lambda ) \ ((J*Phi).' * Y);
% alpha = ( (J*Phi).' * (J*Phi) + gamma1*I + gamma2*Phi.'*L*Phi ) \ ((J*Phi).' * Y);

alpha = (J*(Phi*Lambda*Phi.') + gamma1*I + gamma2*L*(Phi*Lambda*Phi.'))\Y;
b = 0;
end
