function sPreset = UpdatePreset(sPreset, b_interpEigenvecs)
sPreset.nSignals = b_interpEigenvecs*sPreset.M + (1-b_interpEigenvecs)*sPreset.nSignals;
sPreset.b_compareMethods = ~b_interpEigenvecs && sPreset.b_compareMethods;
sPreset.b_interpEigenvecs = b_interpEigenvecs;
sPreset.b_runGraphSignals = ~b_interpEigenvecs && sPreset.b_runGraphSignals;
end