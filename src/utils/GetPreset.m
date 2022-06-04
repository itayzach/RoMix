function sPreset = GetPreset(presetName)
if isstruct(presetName)
    sPreset = presetName;
elseif ischar(presetName)
    sPreset = eval(presetName);
end
PrintPresetParamsLatex(sPreset);
PrintPreset(sPreset);
end