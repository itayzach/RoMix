
    
reg1 = Phi.'*L*Phi;
n = size(Phi,1);
I = eye(n);
lambda = 1./diag(invLambda);
r = diag(exp(-lambda));
reg2 = (Phi.'*Phi)*r*(Phi.'*Phi); %Phi.'*(I - Phi*r*Phi.')*Phi;
matToInv = JPhi.'*JPhi + gamma1*invLambda + gamma2*(reg2);