function [] = VerifyMercerTheorem(sSimParams)
% Verify k(x,y) = sum_{m=0}^inf ~ sum_{m=0}^M( lambda_m*phi_m(x)*phi_m(y) )

if sSimParams.dim == 1
    nyPoints = 1000;
    yMax = 3;
    yMin = -2;
    
    nxPoints = 1000;
    xMax = 3;
    xMin = -2;
    
    y = zeros(nyPoints, sSimParams.dim);
    x = zeros(nxPoints, sSimParams.dim);
    for d = 1:sSimParams.dim
        y(:, d) = (yMax - yMin)*rand(nyPoints, 1) + yMin;
        x(:, d) = (xMax - xMin)*rand(nxPoints, 1) + xMin;
    end
    
    nxPoints = size(x,1);
    nyPoints = size(y,1);
    
else
    if sSimParams.twomoons_dataset
        x = sSimParams.sDataset.x;
        if sSimParams.twomoons_scale
            x = 10*x;
        end
        xMax = max(max(x,[],1));
        xMin = min(min(x,[],1));
        
    else
        nxPoints = 200;
        xMax = 5; %25;
        xMin = -5; %-13;
        x = (xMax - xMin)*rand(nxPoints, 2) + xMin;
    end
    step = (xMax - xMin)/100;
    x1 = xMin:step:xMax;
    x2 = x1;
    [XX1,XX2] = meshgrid(x1,x2);
    
    y = [XX1(:) XX2(:)];
    
    nxPoints = length(x);
    nyPoints = length(y);
end

mPhi_x = zeros(nxPoints, sSimParams.MercerM);
mPhi_y = zeros(nyPoints, sSimParams.MercerM);
vLambda = zeros(1, sSimParams.MercerM);
for i = 0:sSimParams.MercerM-1
    if sSimParams.dim == 1
        m = i;
    elseif sSimParams.dim == 2
        m = OneDim2TwoDimIndex(sSimParams.multindexToSingleIndexMap(i+1)-1);
    else
        error('conversion from i to m d-dim index is not implemented');
    end
    vLambda(i+1) = lambda(sSimParams, m);
    mPhi_x(:, i+1) = phi(sSimParams, m, x);
    mPhi_y(:, i+1) = phi(sSimParams, m, y);
%     if sSimParams.dim == 2
%         if m(2) == 0
%             fprintf('-------------------------------------------------\n');
%         end
%         fprintf('%d\t:\t%.16f\t[%d\t%d]\n', i, vLambda(i+1), m(1), m(2));
%     end
end
fprintf('last eigenvalue is vLambda(%d) = %.12f\n', length(vLambda), vLambda(end));

if sSimParams.constsType == 1
    lhs = calckernel('rbf', sSimParams.ell, y, x);
elseif sSimParams.constsType == 2
    lhs = calckernel('rbf', sSimParams.omega, y, x);
else
    error('');
end
rhs = mPhi_x * diag(vLambda) * mPhi_y.';

figure;
subplot(2,1,1);
    imagesc(lhs); colorbar;
    title('kernel');
    
subplot(2,1,2);
    imagesc(rhs); colorbar;
    title('mercer');


isalmostequal(lhs, rhs, 1e-10, sprintf('Mercer''s theorem failed. Perhaps M = %d is not enough?', sSimParams.MercerM));
fprintf('Mercer''s theorem confirmed for M = %d\n', sSimParams.MercerM);
end