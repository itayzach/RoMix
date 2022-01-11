function [vPhi_m] = phi(sKernelParams, c, m, x)
    assert(length(m) == size(x,2))
    vPhi_m = ones(size(x,1),1); 
    for d = 1:length(m)
        if sKernelParams.sDistParams.sigma{c}(d) == 0
            continue
        end
        vPhi_m = vPhi_m .* phi_d(sKernelParams, c, m(d), x*sKernelParams.sDistParams.u{c}(:,d), d);
    end
end