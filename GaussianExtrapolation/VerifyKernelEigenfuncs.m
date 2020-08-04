function [] = VerifyKernelEigenfuncs(sSimParams)

assert(sSimParams.dim <= 2);
%% Check eigenfunctions

nPoints = 10;
yMax = sSimParams.xMax;
yMin = sSimParams.xMin;

y = zeros(nPoints, sSimParams.dim);
% x = zeros(length(xAxis1d), sSimParams.dim);
for d = 1:sSimParams.dim
    y(:, d) = (yMax - yMin)*rand(nPoints, 1) + yMin;
%     x(:, d) = xAxis1d;
end

for i = 0:sSimParams.RkhsM-1
    if sSimParams.dim == 2
        m = OneDim2TwoDimIndex(sSimParams.multindexToSingleIndexMap(i+1)-1);
    else
        m = i;
    end
    
    %% rhs
    lambda_m = lambda(sSimParams, m);
    if lambda_m < 1e-20
        fprintf('VerifyRKHS: lambda_m < 1e-20, breaking...\n');
        break;
    end
    
    vPhi_m_y = phi(sSimParams, m, y);
    rhs = lambda_m * vPhi_m_y;
    
    %% lhs
    if sSimParams.dim == 1
        lhs = zeros(nPoints,1);
        for j = 1:nPoints
            integrand = @(x) kernel(sSimParams, y(j), x).*phi(sSimParams, m, x).*p(sSimParams, x);
            lhs(j) = integral(integrand,-1e3,1e3,'ArrayValued',true);
        end
    else
        warning('VerifyRKHS does not support multi-index yet...')
        return;
        lhs_d = zeros(nPoints,sSimParams.dim);
        rhs_d = zeros(nPoints,sSimParams.dim);
        for d = 1:sSimParams.dim
            vPhi_m_x = phi_d(sSimParams, i, xAxis1d, d);
            vP_x = p(sSimParams, xAxis1d, d);

            % Make sure that phi_m is an eigenfunction of the kernel by:
            % (for a fixed y)
            %   rhs = lambda_m * phi_m_(y)
            %   lhs = <Ky, phi_m> = integral_x( Ky(x)phi_m(x)p(x)dx )

            vPhi_m_y = phi_d(sSimParams, i, y(:, d), d);
            rhs_d(:, d) = lambda_m * vPhi_m_y;

            for j = 1:nPoints
                vKernel_y_x = kernel(sSimParams, y(j, d), xAxis1d);
                integral_1d = sum(vKernel_y_x.*vPhi_m_x.*vP_x*dx); %1-D integral over single y_{j,d}
                lhs_d(j, d) = integral_1d;
            end
        end
        rhs = prod(rhs_d, 2);
        lhs = prod(lhs_d, 2); % multiply all D integrals
    end
    fprintf('m = %d\n', i)
    fprintf('y            = ')
    fprintf('%10.6f  ', y.');
    fprintf('\n');
    fprintf('lambda * phi = ')
    fprintf('%10.6f  ', rhs.');
    fprintf('\n');
    fprintf('<Ky, phi>    = ')
    fprintf('%10.6f  ', lhs.');
    fprintf('\n');
    isalmostequal(rhs, lhs, 1e-12, sprintf('m = %d failed...', i))
    fprintf('(lhs) <Ky, phi_m> = lambda_m*phi_m(y) (rhs) confirmed\n');
    fprintf('\n');
end

end