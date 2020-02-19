function [] = VerifyMercerTheorem(sParams)
% Verify k(x,y) = sum_{m=0}^inf ~ sum_{m=0}^M( lambda_m*phi_m(x)*phi_m(y) )

nyPoints = 10;
yMax = 2;
yMin = -2;

nxPoints = 10;
xMax = 2;
xMin = -2;

y = zeros(nyPoints, sParams.dim);
x = zeros(nxPoints, sParams.dim);
for d = 1:sParams.dim
    y(:, d) = (yMax - yMin)*rand(nxPoints, 1) + yMin;
    x(:, d) = (xMax - xMin)*rand(nyPoints, 1) + xMin;
end

nxPoints = size(x,1);
nyPoints = size(y,1);

lhs = zeros(nxPoints, nyPoints);
rhs = zeros(nxPoints, nyPoints);

for i = 1:length(x)
    for j = 1:length(y)
        vLambda_Phix_Phiy = zeros(1, sParams.dim);
        for d = 1:sParams.dim
            for m = 0:sParams.M-1    
                [phi_x_d, lambda] = phi(sParams.a, sParams.b, m, x(i,d));
                [phi_y_d, ~] = phi(sParams.a, sParams.b, m, y(j,d));
                
                vLambda_Phix_Phiy(d) = vLambda_Phix_Phiy(d) + lambda*phi_x_d*phi_y_d;
            end
        end
        rhs(i,j) = prod(vLambda_Phix_Phiy, 2);
        k_d = kernel(x(i,:), y(j,:), sParams.l);
        lhs(i,j) = prod(k_d, 2);
    end
end


isalmostequal(rhs, lhs, 1e-10, sprintf('Mercer''s theorem failed. Perhaps M = %d is not enough?', sParams.M));
fprintf('Mercer''s theorem confirmed\n');
end