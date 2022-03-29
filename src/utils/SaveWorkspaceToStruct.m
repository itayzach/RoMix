function sWorkspace = SaveWorkspaceToStruct()
cWorkspace = evalin('caller', 'who');
for ind = 1:numel(cWorkspace)
    curVar = evalin('caller', cWorkspace{ind});
    sWorkspace.(cWorkspace{ind}) = curVar;
end
end