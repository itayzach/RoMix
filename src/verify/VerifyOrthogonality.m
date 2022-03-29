function [] = VerifyOrthogonality(sSimParams)
dx = 0.01;
x_1d_axis = (-1e3:dx:1e3-dx).';

mPhi = zeros(length(x_1d_axis), sSimParams.dim, sSimParams.OrthM);
for m = 0:sSimParams.OrthM-1
    for d = 1:sSimParams.dim
        vPhi_m_x = phi_d(sSimParams, m, x_1d_axis, d);
        mPhi(:, d, m+1) = vPhi_m_x;
    end
end

I = zeros(sSimParams.OrthM);
vP_x = p(sSimParams, x_1d_axis);
for k = 0:sSimParams.OrthM-1
    for m = 0:sSimParams.OrthM-1
        % Following line is 
        %   <phi_m, phi_k> = 
        %   integral_x1( phi_m(x1)phi_k(x1)p(x1)dx1 ) *
        %   integral_x2( phi_m(x2)phi_k(x2)p(x2)dx2 ) * 
        %   * ... *
        %   integral_xd( phi_m(xd)phi_k(xd)p(xd)dxd ) *
        vPhim_Phik_d = sum(mPhi(:, :, m+1).*mPhi(:, :, k+1).*vP_x*dx);
        I(m+1,k+1) = prod(vPhim_Phik_d, 2);
        if abs(I(m+1,k+1)) < 1e-10
            I(m+1,k+1) = 0;
        end
        if abs(I(m+1,k+1) - 1) < 1e-10
            I(m+1,k+1) = 1;
        end
    end
end
isalmostequal(I, eye(sSimParams.OrthM), 0, 'eigenfunctions are not orthonormal');
fprintf('All %d eigenfunctions are orthonormal\n', sSimParams.OrthM);
end