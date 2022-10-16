function PrintAccuracy(sPreset, vAccRec, vAccRecStd, vAccInt, vAccIntStd, vTrainTime, vIntTime, cMethods)
assert(sPreset.nSignals)
nMethods = numel(cMethods);

fprintf('%s\n',repelem('-',15*5+7));
fprintf('%-15s| %-16s| %-16s| %-15s| %-15s\n', 'Method','Interp. Acc.','Extrap. Acc.','Train time [s]','Extrap. time [s]');
fprintf('%s\n',repelem('-',15*5+7));
for ii = 1:nMethods
    if sPreset.R > 0
        fprintf('%-15s| %.4f ± %-7.2f| %.4f ± %-8.2f| %-15.2f| %-15.2f\n',cMethods{ii},vAccRec(ii),vAccRecStd(ii),vAccInt(ii),vAccIntStd(ii),vTrainTime(ii),vIntTime(ii));
    else
        fprintf('%-15s| %6-15.4f| %-15.2f| %-15.2f| %-15.2f\n',cMethods{ii},vAccRec(ii),vAccInt(ii),vTrainTime(ii),vIntTime(ii));
        %s = [s, cMethods{ii} ' : ' num2str(vAccRec(ii), '%.2f') ' | ', num2str(vAccInt(ii), '%.2f'), newline];
    end
    fprintf('%s\n',repelem('-',15*5+7));
end

end

function spacedStr = AddSpaceBeforeUpper(str)
isUpperCase = isstrprop(str, 'upper');
indVec = 1:numel(str);
if sum(isUpperCase) == numel(str) || sum(isUpperCase) == 1
    spacedStr = str;
else
    upperInd = indVec(isUpperCase);
    spacedStr = [str(1:upperInd(2)-1), ' ', str(upperInd(2):end)];
end

end