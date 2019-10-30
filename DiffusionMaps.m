%% Init
clc; clear; close all;
rng('default');

%%
W = [0 4 2 2;
     4 0 8 0;
     2 8 0 0;
     2 0 0 0];
W = [0 1 0 1;
     1 0 1 0;
     0 1 0 1;
     1 0 1 0];
 
X_test = randn(4, 4);

gaussSigma = 1;
W = exp(-(pdist2(X_test, X_test).^2)./(2*gaussSigma^2));
 
D = diag(sum(W));

P = D^-1 * W;
L = D - W;
Ln = (D^-0.5)*L*(D^-0.5);
I_N = eye(length(W));

T = 5;
Pt = zeros(T, length(P), length(P));
Phi_Pt = zeros(T, length(P), length(P));
Lambda_Pt = zeros(T, length(P), length(P));
Psi_Pt = zeros(T, length(P), length(P));


for i = 1:T
    if i == 1
        Pt(1, :, :) = P;
    else
        Pt(i, :, :) = squeeze(Pt(i-1, :, :)) * P;
    end
    % Psi are the right eigenvectors ( P*Psi = Psi*Lambda   ==> P = Psi*Lambda*inv(Psi)    )
    % Lambda are the eigenvalues
    % Phi are the left eigenvectors  ( Phi'*P = Lambda*Phi' ==> P =  inv(Phi')*Lambda*Phi' )
    [Psi_Pt(i, :, :), Lambda_Pt(i, :, :), Phi_Pt(i, :, :)] = eig(squeeze(Pt(i,:,:)));
    [Lambda_P, idx_P] = sort(diag(squeeze(Lambda_Pt(i, :, :))), 'descend');
    Lambda_Pt(i, :, :) = diag(Lambda_P);
    Psi_P = squeeze(Psi_Pt(i, :, :));
    Phi_P = squeeze(Phi_Pt(i, :, :));
    Psi_Pt(i, :, :) = Psi_P(:, idx_P);
    Phi_Pt(i, :, :) = Phi_P(:, idx_P);
end

% [Phi, Lambda, Psi] = svd(P);


Psi_P = squeeze(Psi_Pt(1, :, :));
Lambda_P = squeeze(Lambda_Pt(1, :, :));
Phi_P = squeeze(Phi_Pt(1, :, :));


[Psi_Ln, Lambda_Ln] = eig(Ln);
[Lambda_Ln, idx_Ln] = sort(diag(Lambda_Ln), 'ascend');
Lambda_Ln = diag(Lambda_Ln);
Psi_Ln = Psi_Ln(:, idx_Ln);

[Psi_L, Lambda_L] = eig(L);
[Lambda_L, idx_L] = sort(diag(Lambda_L), 'ascend');
Lambda_L = diag(Lambda_L);
Psi_L = Psi_L(:, idx_L);

assert(all(diag(Lambda_Ln)) >= 0);                                % Ln is PSD, so all eigenvalues must be >= 0
isalmostequal(Psi_Ln*Lambda_Ln*Psi_Ln', Ln);          % Ln is symmetric
isalmostequal(Psi_L*Lambda_L*Psi_L', L);              % L is symmetric
isalmostequal(I_N - Lambda_Ln, Lambda_P);      % 
isalmostequal(Phi_P'^-1*Lambda_P*Phi_P', P);      % Left eigenvectors of P
isalmostequal(Psi_P*Lambda_P*Psi_P^-1, P);        % Right eigenvectors of P
D^-0.5 * Psi_Ln

keyboard