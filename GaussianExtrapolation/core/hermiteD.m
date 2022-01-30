function res = hermiteD(n, beta,xu,mu,sigma)
[N, D] = size(xu);

if D > 200
    zerosInd = find(n == 0);
    onesInd = find(n == 1);
    twosInd = find(n == 2);
    res = ones(N,D);
    if size(zerosInd,2) == D
        return
    end
    x_onesInd = hermiteArg(beta(onesInd),xu(:,onesInd),mu(onesInd),sigma(onesInd));
    res(:,onesInd) = 2*x_onesInd;
    if ~isempty(twosInd)
        x_twosInd = hermiteArg(beta(twosInd),xu(:,twosInd),mu(twosInd),sigma(twosInd));
        res(:,twosInd) = 4*x_twosInd(twosInd).^2 - 2;
    end
    assert(size([zerosInd, onesInd, twosInd],2) == D)
else
    P = max(n)+1;
    p = P-1:-1:0;
    % p = 0:P-1;
    %% Get Hermite coefficients
    h = LoadHermiteCoeffs(n);
    %% Calculate polynomial
    x = hermiteArg(beta,xu,mu,sigma);
    X = permute(repmat(x,[1, 1, P]), [1 3 2]);
    H = reshape(h, [1, P, D]);
    Xp_h = pagemtimes(X.^p,'none',H,'transpose');
    res = reshape(Xp_h, N, D);

end


end

function h = LoadHermiteCoeffs(n)
% h should be PxD
P = max(n)+1;
D = size(n,2);
n_fact = factorial(n);
h = zeros(P,D);
for d=1:D
    m = 0:floor(n(d)/2);
    rem = mod(n(d),2);
    switch n(d)
        case 0
            coeffs = 1; % 1*x^0
        case 1
            coeffs = 2; % 2*x^1
        case 2
            coeffs = [4, -2]; % 4*x^2 -2*x^0
        case 3
            coeffs = [8, -12];
        otherwise
            coeffs = n_fact(d) .* (-1).^m ./ (factorial(m) .* factorial(n(d)-2.*m)) .* 2.^(n(d)-2.*m);
    end
    h(2*m+1+rem,d) = fliplr(coeffs);
end
h = flipud(h);
end