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
        xMax = 1; %25;
        xMin = -1; %-13;
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

tPhi_x = zeros(nxPoints, sParams.MercerM, sParams.dim);
tPhi_y = zeros(nyPoints, sParams.MercerM, sParams.dim);
vLambda = zeros(1, sParams.MercerM, sParams.dim);
for m = 0:sParams.MercerM-1
    vLambda(1, m+1, :) = lambda(sParams, m);
%     if vLambda(m+1) < 10e-16
%         fprintf('VerifyMercerTheorem: lambda_m < 1e-20, breaking...\n');
%         break
%     end
    assert(~any(isnan(vLambda(m+1))));
    for d = 1:sParams.dim
        tPhi_x(:,m+1,d) = phi(sParams, m, x(:,d), d);
        assert(~any(isnan(squeeze(tPhi_x(:,m+1,d)))));
        tPhi_y(:,m+1,d) = phi(sParams, m, y(:,d), d);
        assert(~any(isnan(squeeze(tPhi_y(:,m+1,d)))));
    end
end




for i = 1:length(x)
    for j = 1:length(y)
        vLambda_Phix_Phiy = sum(vLambda(1,:,:).*tPhi_x(i,:,:).*tPhi_y(j,:,:), 2);
        rhs(i,j) = prod(vLambda_Phix_Phiy, 3);
        k_d = kernel(sParams, x(i,:), y(j,:));
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

isalmostequal(lhs, rhs, 1e-10, sprintf('Mercer''s theorem failed. Perhaps M = %d is not enough?', sParams.MercerM));
fprintf('Mercer''s theorem confirmed for M = %d\n', sParams.MercerM);
end