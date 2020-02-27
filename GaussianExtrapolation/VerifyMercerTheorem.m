function [] = VerifyMercerTheorem(sParams, sSimParams)
% Verify k(x,y) = sum_{m=0}^inf ~ sum_{m=0}^M( lambda_m*phi_m(x)*phi_m(y) )

if sParams.dim == 1
    nyPoints = 1000;
    yMax = 3;
    yMin = -2;
    
    nxPoints = 1000;
    xMax = 3;
    xMin = -2;
    
    y = zeros(nyPoints, sParams.dim);
    x = zeros(nxPoints, sParams.dim);
    for d = 1:sParams.dim
        y(:, d) = (yMax - yMin)*rand(nyPoints, 1) + yMin;
        x(:, d) = (xMax - xMin)*rand(nxPoints, 1) + xMin;
    end
    
    nxPoints = size(x,1);
    nyPoints = size(y,1);
    
else
    if sSimParams.twomoons_dataset
        twomoons = load('2moons.mat', 'x');
        x = twomoons.x;
        if sSimParams.twomoons_scale
            x = 10*x;
        end
        xMax = max(max(x,[],1));
        xMin = min(min(x,[],1));
        
    else
        nxPoints = 200;
        if sSimParams.twomoons_scale
            xMax = 5;25;
            xMin = -5;-13;
        else
            xMax = 3;
            xMin = -2;
        end
        x = (xMax - xMin)*rand(nxPoints, 2) + xMin;
    end
    step = (xMax - xMin)/100;
    x1 = xMin:step:xMax;
    x2 = x1;
    [XX1,XX2] = meshgrid(x1,x2);
    
    y = [XX1(:) XX2(:)];
    
    nxPoints = length(x);
    nyPoints = length(y);
    lhs = zeros(nxPoints, nyPoints);
    rhs = zeros(nxPoints, nyPoints);
end

tPhi_x = zeros(nxPoints, sParams.M, sParams.dim);
tPhi_y = zeros(nyPoints, sParams.M, sParams.dim);
vLambda = zeros(1, sParams.M);
for m = 0:sParams.M-1
    for d = 1:sParams.dim
        [tPhi_x(:,m+1,d), vLambda(m+1)] = phi(sParams.a, sParams.b, m, x(:,d));
        [tPhi_y(:,m+1,d), ~] = phi(sParams.a, sParams.b, m, y(:,d));
    end
end




for i = 1:length(x)
    for j = 1:length(y)
        vLambda_Phix_Phiy = sum(vLambda.*tPhi_x(i,:,:).*tPhi_y(j,:,:), 2);
        rhs(i,j) = prod(vLambda_Phix_Phiy, 3);
        k_d = kernel(x(i,:), y(j,:), sParams.l);
        lhs(i,j) = prod(k_d, 2);
    end
end

if sParams.dim == 2
    figure;
    subplot(2,1,1);
    imagesc(lhs); colorbar;
    title('kernel');
    
    subplot(2,1,2);
    imagesc(rhs); colorbar;
    title('mercer');
end

isalmostequal(rhs, lhs, 1e-10, sprintf('Mercer''s theorem failed. Perhaps M = %d is not enough?', sParams.M));
fprintf('Mercer''s theorem confirmed\n');
end