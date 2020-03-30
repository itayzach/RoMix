[~, scriptName, ~] = fileparts(mfilename('fullpath'));
fprintf('[%s] Adding paths...\n', scriptName);
addpath(genpath('GaussianExtrapolation'));
addpath(genpath('ManifoldLearn'));
addpath('data');
fprintf('[%s] Ready.\n', scriptName);