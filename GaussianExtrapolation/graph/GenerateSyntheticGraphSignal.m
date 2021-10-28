function [sig, sigRef] = GenerateSyntheticGraphSignal(V, VRef)

[n, M] = size(V);
N = length(VRef);
interpRatio = N/n;

gamma = 0.001*M;
sigHat = exp(-gamma*(1:M)');
sig = 100*V*sigHat;
sigRef = 100*sqrt(interpRatio)*VRef*sigHat;

% noise = 0*randn(n,1);
% sigNoisy = sig + noise(1:n);
% sigRef = sqrt(interpRatio)*VRef*sigHat;
% sigRef = VRef*sigHat;
% sigRefNoisy = sigRef + noise;
end
