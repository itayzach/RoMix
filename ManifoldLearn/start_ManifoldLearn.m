fprintf('Starting ManifoldLearn...\n')
addpath('data');
if isempty(dir('*.mex*'))
    fprintf('Compiling mex files for SVM...\n')
    mex -O -c svmprecomputed.cpp
    mex -O mexGramSVMTrain.cpp  svmprecomputed.obj
end
fprintf('ManifoldLearn is ready\n')