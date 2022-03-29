function ReadWorkspaceStructToBase(sWorkspace)
    fieldNames = fieldnames(sWorkspace);
    for i = 1:numel(fieldNames)
        varName = string(fieldNames(i));
        varValue = sWorkspace.(varName);
%         if (isstruct(field))
%             unpackStruct(field);
%             continue;
%         end
        assignin('base', varName, varValue);
    end
end