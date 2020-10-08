function [f,f_hat] = GenerateGraphSignal(G, mPhiAnalytic)

% this:
% paramf.log = 1;
% Nf = 5;
% g = gsp_design_warped_translates(G, Nf,paramf);  
% s = sign(G.U(:,2));
% sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
% f = sf*[0 1 0 0 0]';
% f_hat = G.U.'*f;

% or this:
% orig_c = [0 0 100 0 0 0 0 0 0 0]';
% f = G.U(:,1:nEigs)*orig_c;

% or this:
% K = floor(samplingRatio*G.N);
% f_hat_K = 100*exp(-0.5*(1:K))';
% f_hat = [f_hat_K; zeros(G.N-K,1)];
% f = G.U*(f_hat);

% or this:
% f_hat = 5*exp(-0.5*(1:G.N))';
% f = G.U*(f_hat);

% or this:
k0 = round(0.01*G.N);
f_hat = zeros(G.N,1);
f_hat(1:k0) = 5*sort(abs(randn(k0,1)), 'descend');
f = G.U*f_hat;
% f = V.'*f_hat;

% or this:
% nEigs = size(mPhiAnalytic,2);
% c_orig = randn(nEigs, 1);
% f = mPhiAnalytic*c_orig;
% f_hat = G.U*f;


end

