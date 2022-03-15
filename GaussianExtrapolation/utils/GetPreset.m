function sPreset = GetPreset(presetName)
if isstruct(presetName)
    sPreset = presetName;
elseif ischar(presetName)
    sPreset = eval(presetName);
end
PrintPreset(sPreset);
end