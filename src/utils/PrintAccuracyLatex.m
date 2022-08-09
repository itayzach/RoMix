function PrintAccuracyLatex(sPreset, vAccRec, vAccRecStd, vAccInt, vAccIntStd, vTrainTime, cMethods)
assert(sPreset.nSignals == 1)
nMethods = numel(cMethods);
s = [];
%s = [s, '\begin{table}[htbp]', newline ];
s = [s, '\begin{minipage}[c]{0.5\textwidth}', newline ];
s = [s, '    \centering', newline ];
s = [s, '    \begin{tabular}{ |l|l|l|l| }', newline ];
s = [s, '    \hline', newline ];
s = [s, '    Method & Interp. Acc. & Extrap. Acc. & Train Time [s] \\', newline ];
s = [s, '    \hline\hline', newline ];
for ii = 1:nMethods
    if sPreset.R > 1
        s = [s, '    ', cMethods{ii} ' &  $' num2str(vAccRec(ii), '%.2f') ' \pm ' num2str(vAccRecStd(ii), '%.2f') '\% $', ...
            ' &  $' num2str(vAccInt(ii), '%.2f') ' \pm ' num2str(vAccIntStd(ii), '%.2f') ' \%$ & $ ', num2str(vTrainTime(ii), '%.2f'), ' $ \\', newline ];
    else
        s = [s, '    ', cMethods{ii} ' &  $', num2str(vAccRec(ii), '%.2f'), '\% $', ...
            ' &  $', num2str(vAccInt(ii), '%.2f'), ' \%$ & $ ', num2str(vTrainTime(ii), '%.2f'),  ' $ \\', newline];
    end
    s = [s, '    \hline', newline ];
end
s = [s, '    \end{tabular}', newline ];
%s = [s, '    \captionof{table}{Accuracy for ' AddSpaceBeforeUpper(sPreset.verticesPDF) ' dataset}', newline ];
%s = [s, '    \label{tab:acc_' sPreset.verticesPDF '}', newline ];
s = [s, '\end{minipage}'];
%s = [s, '\end{table}'];
fprintf('********************** COPY BEGIN **********************\n')
fprintf('%s\n', s);
fprintf('*********************** COPY END ***********************\n')
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