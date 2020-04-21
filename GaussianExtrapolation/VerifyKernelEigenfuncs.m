function [] = VerifyKernelEigenfuncs(sParams)

assert(sParams.dim <= 2);
%% Check eigenfunctions
dx = 0.05;
xAxis1d = (-1e3:dx:1e3-dx).';
nPoints = 10;
yMax = 3;
yMin = -3;

y = zeros(nPoints, sParams.dim);
% x = zeros(length(xAxis1d), sParams.dim);
for d = 1:sParams.dim
    y(:, d) = (yMax - yMin)*rand(nPoints, 1) + yMin;
%     x(:, d) = xAxis1d;
end

for i = 0:sParams.RkhsM-1
    if sParams.dim == 2
        m = OneDim2TwoDimIndex(i);
    else
        m = i;
    end
    
    %% rhs
    lambda_m = lambda(sParams, m);
    if lambda_m < 1e-20
        fprintf('VerifyRKHS: lambda_m < 1e-20, breaking...\n');
        break;
    end
    
    vPhi_m_y = phi(sParams, m, y);
    rhs = lambda_m * vPhi_m_y;
    
    %% lhs
    if sParams.dim == 1
        lhs = zeros(nPoints,1);
        for j = 1:nPoints
            integrand = @(x) kernel(sParams, y(j), x).*phi(sParams, m, x).*p(sParams, x);
            lhs(j) = integral(integrand,-1e3,1e3,'ArrayValued',true);
        end
    else
        warning('VerifyRKHS does not support multi-index yet...')
        return;
        lhs_d = zeros(nPoints,sParams.dim);
        rhs_d = zeros(nPoints,sParams.dim);
        for d = 1:sParams.dim
            vPhi_m_x = phi_d(sParams, i, xAxis1d, d);
            vP_x = p(sParams, xAxis1d, d);

            % Make sure that phi_m is an eigenfunction of the kernel by:
            % (for a fixed y)
            %   rhs = lambda_m * phi_m_(y)
            %   lhs = <Ky, phi_m> = integral_x( Ky(x)phi_m(x)p(x)dx )

            vPhi_m_y = phi_d(sParams, i, y(:, d), d);
            rhs_d(:, d) = lambda_m * vPhi_m_y;

            for j = 1:nPoints
                vKernel_y_x = kernel(sParams, y(j, d), xAxis1d);
                integral_1d = sum(vKernel_y_x.*vPhi_m_x.*vP_x*dx); %1-D integral over single y_{j,d}
                lhs_d(j, d) = integral_1d;
            end
        end
        rhs = prod(rhs_d, 2);
        lhs = prod(lhs_d, 2); % multiply all D integrals
    end
    fprintf('m = %d\n', i)
    fprintf('lambda * phi = ')
    fprintf('%.6f  ', rhs.');
    fprintf('\n');
    fprintf('<Ky, phi>     = ')
    fprintf('%.6f  ', lhs.');
    fprintf('\n');
    isalmostequal(rhs, lhs, 1e-12, sprintf('m = %d failed...', i))
    fprintf('(lhs) <Ky, phi_m> = lambda_m*phi_m(y) (rhs) confirmed\n');
    fprintf('\n');
end

end