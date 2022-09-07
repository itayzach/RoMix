function sPreset = UpdatePreset(sPreset, b_interpEigenvecs, nLabeled)
sPreset.nSignals = b_interpEigenvecs*sPreset.M + (1-b_interpEigenvecs)*sPreset.nSignals;
sPreset.b_compareMethods = ~b_interpEigenvecs && sPreset.b_compareMethods;
sPreset.b_interpEigenvecs = b_interpEigenvecs;
sPreset.b_runGraphSignals = ~b_interpEigenvecs && sPreset.b_runGraphSignals;
sPreset.nLabeled = nLabeled;

if b_interpEigenvecs && ~sPreset.b_debugUseAnalytic
    % Numeric eigs() returns the eigenvectors with norm = 1. 
    % When n ~= N, we should renormalize
    sPreset.interpRatio = sPreset.N/sPreset.n;
else
    sPreset.interpRatio = 1;
end

end