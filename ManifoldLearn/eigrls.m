function [alpha, c]=eigrls(Y, Phi, Lambda, gamma1, gamma2, L)

n = size(Y,1);   % total examples
l = sum(abs(Y)); % labeled examples
u = n - l;       % unlabeled examples

M = size(Phi,2);
I_n = eye(n);
J = diag(abs(Y));


% alpha = (J*(Phi*Lambda*Phi.') + gamma1*I_n + gamma2*L*(Phi*Lambda*Phi.'))\Y;
alpha = 0;

invLambda = diag(diag(1./Lambda));
c = ( (J*Phi).' * (J*Phi) + gamma1*invLambda + gamma2*Phi.'*L*Phi ) \ ((J*Phi).' * Y);


end
