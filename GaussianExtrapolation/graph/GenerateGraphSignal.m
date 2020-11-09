function [f, coeffs] = GenerateGraphSignal(B, graphSignalModel)

N = size(B,1);

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
    
    k0 = round(0.01*N);
    f_hat = zeros(N,1);
    f_hat(1:k0) = 5*sort(abs(randn(k0,1)), 'descend');
    f = B*f_hat;
    coeffs = f_hat;
else
    if strcmpi(graphSignalModel, 'V_c') || strcmpi(graphSignalModel, 'Phi_c')
        nEigs = size(B,2);
        c_orig = 10*randn(nEigs, 1);
        f = B*c_orig;
        coeffs = c_orig;
    elseif strcmpi(graphSignalModel, 'K_alpha')
        alpha = randn(N,1);
        f = B*alpha; % K*alpha
        coeffs = alpha;
    else
        error('unknown graphSignalModel');
    end
end


end

