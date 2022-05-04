function PrintPreset(sPreset)
d = sPreset.dim;
k = sPreset.gmmNumComponents;
fprintf('*********************************************************\n')
fprintf('* Running: %d-D %s\n', sPreset.dim, sPreset.verticesPDF)
fprintf('*********************************************************\n')
fprintf('* Number of GMM parameters : \n\t(dxd + d)*k = (%dx%d + %d)*%d = %d\n', ...
    d,d,d,k,(d*d + d)*k)
fprintf('* Number of coefficients : %d\n', sPreset.MTilde)
fprintf('*********************************************************\n')
end