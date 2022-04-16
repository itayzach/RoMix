function PrintPreset(sPreset)
d = sPreset.dim;
k = sPreset.gmmNumComponents;
fprintf('*********************************************************\n')
fprintf('* Running: %d-D %s\n', sPreset.dim, sPreset.verticesPDF)
fprintf('*********************************************************\n')
fprintf('* Number of parameters (mu,sigma,k): \n\t(dxd + d)*k = (%dx%d + %d)*%d = %d\n', ...
    d,d,d,k,(d*d + d)*k)
fprintf('*********************************************************\n')
end