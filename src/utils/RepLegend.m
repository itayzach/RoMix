function cellArray = RepLegend(strName, plotInd)
% Run example: RepLegend('v^{{\\bf ref}}', [2, 3, 4]) 
%     ---> {'$v^{{\bf ref}}_2$'}    {'$v^{{\bf ref}}_3$'}    {'$v^{{\bf ref}}_4$'}
cellArray = strsplit(sprintf(['$', strName, '_{%d}$,'], plotInd), ',');
cellArray = cellArray(~cellfun(@isempty, cellArray));
end