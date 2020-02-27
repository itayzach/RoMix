function [] = VerifyRHKS(sParams)
%% Check eigenfunctions
dx = 0.01;
x = (-1e3:dx:1e3-dx).';

nPoints = 10;
yMax = 1;
yMin = -1;

y = zeros(nPoints, sParams.dim);
for d = 1:sParams.dim
    y(:, d) = (yMax - yMin)*rand(nPoints, 1) + yMin;
end
rhs = zeros(size(y));
lhs = zeros(size(y));

for m = 0:sParams.M-1
    lambda_m = lambda(sParams.a, sParams.b, m);
    if lambda_m < 1e-20
        fprintf('VerifyRKHS: lambda_m < 1e-20, breaking...\n');
        break;
    end
    for d = 1:sParams.dim
        vPhi_m_x = phi(sParams.a, sParams.b, m, x);
        vP_x = p(x, sParams.sigma);

        % Make sure that phi_m is an eigenfunction of the kernel by:
        % (for a fixed y)
        %   rhs = lambda_m * phi_m_(y)
        %   lhs = <Ky, phi_m> = integral_x( Ky(x)phi_m(x)p(x)dx )
        
        phi_m_y = phi(sParams.a, sParams.b, m, y(:, d));
        rhs(:, d) = lambda_m * phi_m_y;

        for j = 1:nPoints
            vKernel_y_x = kernel(y(j, d), x, sParams.l);
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
    isalmostequal(rhs_d, lhs_d, 1e-12)
    fprintf('(lhs) <Ky, phi_m> = lambda_m*phi_m(y) (rhs) confirmed\n');
    fprintf('\n');
end

end