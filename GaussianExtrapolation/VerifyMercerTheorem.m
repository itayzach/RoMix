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
        mPhi_x = zeros(sParams.M, sParams.dim);
        mPhi_y = zeros(sParams.M, sParams.dim);
        vLambda = zeros(sParams.M, 1);
        for m = 0:sParams.M-1    
            for d = 1:sParams.dim
                [mPhi_x(m+1,d), vLambda(m+1)] = phi(sParams.a, sParams.b, m, x(i,d));
                [mPhi_y(m+1,d), ~] = phi(sParams.a, sParams.b, m, y(j,d));
            end
            
%             vLambda(m+1)*prod(mPhi_x(m+1,:), 2)*prod(mPhi_y(m+1,:), 2)
            rhs(i,j) = rhs(i,j) + vLambda(m+1)*prod(mPhi_x(m+1,:), 2)*prod(mPhi_y(m+1,:), 2);
        end
        
        k_d = kernel(x(i,:), y(j,:), sParams.l);
%         for d = 1:sParams.dim
%             k_d(d) = kernel(x(i,d), y(j,d), sParams.l);
%         end
        lhs(i,j) = prod(k_d, 2);
    end
end


isalmostequal(rhs, lhs, 1e-10);
fprintf('Mercer''s theorem confirmed\n');
end