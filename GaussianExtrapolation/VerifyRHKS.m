function [] = VerifyRHKS(sParams)

assert(sParams.dim <= 2);
%% Check eigenfunctions
dx = 0.005;
x_1d_axis = (-1e3:dx:1e3-dx).';

nPoints = 10;
yMax = 3;
yMin = -3;

y = zeros(nPoints, sParams.dim);
for d = 1:sParams.dim
    y(:, d) = (yMax - yMin)*rand(nPoints, 1) + yMin;
end
rhs = zeros(size(y));
lhs = zeros(size(y));

for i = 0:sParams.RkhsM-1
    if sParams.dim == 2
        m = OneDim2TwoDimIndex(i, sParams.dim);
    else
        m = i;
    end
    
    lambda_m = lambda(sParams, m);
    if lambda_m < 1e-20
        fprintf('VerifyRKHS: lambda_m < 1e-20, breaking...\n');
        break;
    end
    
%     vPhi_m_y = phi(sParams, m, y);
%     rhs_d = lambda_m * vPhi_m_y;
    
    for d = 1:sParams.dim
        vPhi_m_x = phi_d(sParams, i, x_1d_axis, d);
        vP_x = p(sParams, x_1d_axis, d);

        % Make sure that phi_m is an eigenfunction of the kernel by:
        % (for a fixed y)
        %   rhs = lambda_m * phi_m_(y)
        %   lhs = <Ky, phi_m> = integral_x( Ky(x)phi_m(x)p(x)dx )
        
        vPhi_m_y = phi_d(sParams, i, y(:, d), d);
        rhs(:, d) = lambda_m * vPhi_m_y;

        for j = 1:nPoints
            vKernel_y_x = kernel(sParams, y(j, d), x_1d_axis);
            integral_1d = sum(vKernel_y_x.*vPhi_m_x.*vP_x*dx); %1-D integral over single y_{j,d}
            lhs(j, d) = integral_1d;
        end
    end
    rhs_d = prod(rhs, 2);
    lhs_d = prod(lhs, 2); % multiply all D integrals
    fprintf('m = %d\n', i)
    fprintf('lambda * phi = ')
    fprintf('%.6f  ', rhs_d.');
    fprintf('\n');
    fprintf('<Ky, phi>     = ')
    fprintf('%.6f  ', lhs_d.');
    fprintf('\n');
    isalmostequal(rhs_d, lhs_d, 1e-12, sprintf('m = %d failed...', i))
    fprintf('(lhs) <Ky, phi_m> = lambda_m*phi_m(y) (rhs) confirmed\n');
    fprintf('\n');
end

end