function PrintAccuracy(sPreset, vAccRec, vAccRecStd, vAccInt, vAccIntStd, vTrainTime, vIntTime, cMethods)
assert(sPreset.nSignals)
nMethods = numel(cMethods);
methods = [];
accRec = [];
accInt = [];
for ii = 1:nMethods
    methods = [methods, '& ' cMethods{ii} ' ' ];
    if sPreset.R > 0
        accRec = [ accRec, ' &  $' num2str(vAccRec(ii), '%.2f') ' \pm ' num2str(vAccRecStd(ii), '%.2f') '\% $'];
    else
        accRec = [ accRec, ' &  $' num2str(vAccRec(ii), '%.2f') '\% $'];
    end
    if sPreset.R > 0
        accInt = [ accInt, ' &  $' num2str(vAccInt(ii), '%.2f') ' \pm ' num2str(vAccIntStd(ii), '%.2f') ' \%$'];
    else
        accInt = [ accInt, ' &  $' num2str(vAccInt(ii), '%.2f') ' \%$'];
    end
end
fprintf('%s\n',repelem('-',15*5+7));
fprintf('%-15s| %-16s| %-16s| %-15s| %-15s\n', 'Method','Proj. Acc.','Interp. Acc.','Train time','Interp. time');
fprintf('%s\n',repelem('-',15*5+7));
for ii = 1:nMethods
    if sPreset.R > 0
        fprintf('%-15s| %6.2f ± %-7.2f| %.2f ± %-8.2f| %-15.2f| %-15.2f\n',cMethods{ii},vAccRec(ii),vAccIntStd(ii),vAccInt(ii),vAccIntStd(ii),vTrainTime(ii),vIntTime(ii));
    else
        fprintf('%-15s| %6-15.2f| %-15.2f| %-15.2f| %-15.2f\n',cMethods{ii},vAccRec(ii),vAccInt(ii),vTrainTime(ii),vIntTime(ii));
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