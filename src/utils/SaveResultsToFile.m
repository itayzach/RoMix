function SaveResultsToFile(sPreset, sPlotParams, sResults, b_interpEigenvecs)
t = datetime();
if b_interpEigenvecs
    eigOrSig = 'Eigs';
else
    eigOrSig = 'Sig';
end
filename = [num2str(t.Year), num2str(t.Month,'%02d'), num2str(t.Day,'%02d'), '_', ...
    num2str(t.Hour,'%02d'), num2str(t.Minute,'%02d'), num2str(round(t.Second),'%02d'), '_', ...
    sPreset.verticesPDF, '_', eigOrSig, '.mat'];
outputFolder = fullfile('runmat');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
end
save(fullfile(outputFolder,filename), 'sPreset', 'sPlotParams', 'sResults');
end