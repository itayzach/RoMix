function [sig, sigRef] = GenerateSyntheticGraphSignal(V)

[N, M] = size(V);

gamma = 0.001*M;
sigHat = exp(-gamma*(1:M)');
sig = 100*V*sigHat;
sigRef = 100*VRefToCompare*sigHat;

noise = 0*randn(N,1);
sigNoisy = sig + noise(1:n);
sigRef = sqrt(interpRatio)*VRefToCompare*sigHat;
sigRefNoisy = sigRef + noise;
end