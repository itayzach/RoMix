function [f] = GenerateGraphSignal(graphSignalModel, B, coeffs)

if strcmpi(graphSignalModel, 'bandlimited') || strcmpi(graphSignalModel, 'U_fhat')
    % paramf.log = 1;
    % Nf = 5;
    % g = gsp_design_warped_translates(G, Nf,paramf);  
    % s = sign(G.U(:,2));
    % sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
    % f = sf*[0 1 0 0 0]';
    % f_hat = G.U.'*f;

    % orig_c = [0 0 100 0 0 0 0 0 0 0]';
    % f = G.U(:,1:nEigs)*orig_c;

    % K = floor(samplingRatio*N);
    % f_hat_K = 100*exp(-0.5*(1:K))';
    % f_hat = [f_hat_K; zeros(N-K,1)];
    % f = G.U*(f_hat);

    % f_hat = 5*exp(-0.5*(1:N))';
    % f = G.U*(f_hat);

    f = B*coeffs;
elseif strcmpi(graphSignalModel, 'V_c') || strcmpi(graphSignalModel, 'Phi_c')   
    f = B*coeffs;
elseif strcmpi(graphSignalModel, 'K_alpha')
    f = B*coeffs; % K*alpha
else
    error('unknown graphSignalModel');
end


end

