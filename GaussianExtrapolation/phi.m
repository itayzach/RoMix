%% phi (Squared Exponentional)
function [vPhi_m] = phi(sParams, m, x)
    assert(length(m) == size(x,2))
    vPhi_m = ones(size(x,1),1); 
    for d = 1:sParams.dim
        vPhi_m = vPhi_m .* phi_d(sParams, m(d), x(:,d), d);
    end
end