function c = eigrls(Y, Phi, Lambda, gamma1, gamma2, L)

n = size(Y,1);   % total examples
l = sum(abs(Y)); % labeled examples
u = n - l;       % unlabeled examples

M = size(Phi,2);
I_n = eye(n);
J = diag(sign(abs(Y)));
invLambda = diag(diag(1./Lambda));

% Least squares: pinv(A'*A)*(A'*y);
c = ( (J*Phi).' * (J*Phi) + gamma1*invLambda + gamma2*Phi.'*L*Phi ) \ ((J*Phi).' * Y);
% c = pinv((J*Phi)'*(J*Phi)+ gamma1*invLambda + gamma2*Phi.'*L*Phi)*((J*Phi)'*Y);
end
