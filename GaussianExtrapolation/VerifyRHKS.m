function [] = VerifyRHKS(sParams)
%% Check eigenfunctions
dx = 0.005;
x = (-1e3:dx:1e3-dx).';

nPoints = 10;
yMax = 3;
yMin = -3;

y = zeros(nPoints, sParams.dim);
for d = 1:sParams.dim
    y(:, d) = (yMax - yMin)*rand(nPoints, 1) + yMin;
end
rhs = zeros(size(y));
lhs = zeros(size(y));

for m = 0:sParams.RkhsM-1
    lambda_m = lambda(sParams, m);
    if lambda_m < 1e-20
        fprintf('VerifyRKHS: lambda_m < 1e-20, breaking...\n');
        break;
    end
    for d = 1:sParams.dim
        vPhi_m_x = phi(sParams, m, x, d);
        vP_x = p(sParams, x, d);

        % Make sure that phi_m is an eigenfunction of the kernel by:
        % (for a fixed y)
        %   rhs = lambda_m * phi_m_(y)
        %   lhs = <Ky, phi_m> = integral_x( Ky(x)phi_m(x)p(x)dx )
        
        phi_m_y = phi(sParams, m, y(:, d), d);
        rhs(:, d) = lambda_m(d) * phi_m_y;

        for j = 1:nPoints
            vKernel_y_x = kernel(sParams, y(j, d), x);
            integral_1d = sum(vKernel_y_x.*vPhi_m_x.*vP_x*dx); %1-D integral over single y_{j,d}
            lhs(j, d) = integral_1d;
        end
    end
    rhs_d = prod(rhs, 2);
    lhs_d = prod(lhs, 2); % multiply all D integrals
    fprintf('m = %d\n', m)
    fprintf('lambda * phi = ')
    fprintf('%.6f  ', rhs_d.');
    fprintf('\n');
    fprintf('<Ky, phi>     = ')
    fprintf('%.6f  ', lhs_d.');
    fprintf('\n');
    isalmostequal(rhs_d, lhs_d, 1e-12, sprintf('m = %d failed...', m))
    fprintf('(lhs) <Ky, phi_m> = lambda_m*phi_m(y) (rhs) confirmed\n');
    fprintf('\n');
end

end