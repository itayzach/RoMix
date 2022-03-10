function ClearPrevOutputs(outputFolder)
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
else
    % Remove all epsc files
    filePattern = fullfile(outputFolder, '*.eps');
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile(outputFolder, baseFileName);
      delete(fullFileName);
    end
end