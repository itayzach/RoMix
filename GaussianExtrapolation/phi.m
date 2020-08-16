function [vPhi_m] = phi(sKernelParams, m, x)
    assert(length(m) == size(x,2))
    vPhi_m = ones(size(x,1),1); 
    for d = 1:length(m)
        vPhi_m = vPhi_m .* phi_d(sKernelParams, m(d), x*sKernelParams.sDistParams.u(:,d), d);
    end
end