function sig = GenerateSyntheticGraphSignal(x)

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

L = max(x) - min(x);
f = 5/(L);

sig = sin(2*pi*f*x).*exp(-5*L*abs(x-0.5*L));
% sigRef = sin(3*xInt).*exp(-0.5*abs(xInt));

end
