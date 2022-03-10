function [sig, sigRef] = GenerateSyntheticGraphSignal(x, xInt)

% [n, M] = size(V);
% N = length(VRef);
% interpRatio = N/n;
% 
% gamma = 0.001*M;
% sigHat = exp(-gamma*(1:M)');
% sig = 100*V*sigHat;
% sigRef = 100*sqrt(interpRatio)*VRef*sigHat;

% noise = 0*randn(n,1);
% sigNoisy = sig + noise(1:n);
% sigRef = sqrt(interpRatio)*VRef*sigHat;
% sigRef = VRef*sigHat;
% sigRefNoisy = sigRef + noise;

[n,d] = size(xInt);
if d == 1
    L = max(xInt) - min(xInt);
    f = 5/(L);
    sig = sin(2*pi*f.*x).*exp(-5*L.*abs(x-0.5*L));
    sigRef = sin(2*pi*f.*xInt).*exp(-5*L.*abs(xInt-0.5*L));
elseif d == 2
    dims = [1, 2];
    Lx = max(xInt(:, dims)) - min(xInt(:, dims));
    fx = 2./(Lx);
    sig = prod(sin(2*pi*fx.*x(:, dims)),2);
    sigRef = prod(sin(2*pi*fx.*xInt(:, dims)),2);
elseif d == 3
    dims = [1, 3];
    Lx = max(xInt(:, dims)) - min(xInt(:, dims));
    fx = 2./(Lx);
    sig = prod(sin(2*pi*fx.*x(:, dims)),2);
    sigRef = prod(sin(2*pi*fx.*xInt(:, dims)),2);
end
% sigRef = sin(3*xInt).*exp(-0.5*abs(xInt));

end
