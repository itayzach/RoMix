[~, scriptName, ~] = fileparts(mfilename('fullpath'));
fprintf('[%s] Starting...\n', scriptName);
cd ManifoldLearn
if isempty(dir('*.mex*'))
    fprintf('[%s] Compiling mex files for SVM...\n', scriptName)
    mex -O -c svmprecomputed.cpp
    mex -O mexGramSVMTrain.cpp  svmprecomputed.obj
end
fprintf('[%s] Ready.\n', scriptName);
clear scriptName
cd ..