function sPreset = GetPreset(presetName)
if isstruct(presetName)
    sPreset = presetName;
elseif ischar(presetName)
    sPreset = eval(presetName);
end
if ~isfield(sPreset,'nLabeled')
    sPreset.nLabeled = sPreset.n;
end
if ~isfield(sPreset,'cMethods')
    sPreset.cMethods = {'RoMix'};
end
PrintPresetParamsLatex(sPreset);
PrintPreset(sPreset);
end