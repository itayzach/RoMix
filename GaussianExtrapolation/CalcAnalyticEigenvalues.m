function [vLambdaAnaytic, vMultindexToSingleIndexMap] = CalcAnalyticEigenvalues(sSimParams, sKernelParams, dim)

%% Get correct order of eigenvalues (for 1D indexing from multindexing)
if dim == 2
    vLambda_K_before_sort = zeros(max(sSimParams.PlotSpectM,sSimParams.PlotEigenFuncsM),1);
    for i = 0:max(sSimParams.PlotSpectM,sSimParams.PlotEigenFuncsM)-1
        m = OneDim2TwoDimIndex(i);
        vLambda_K_before_sort(i+1) = lambda(sKernelParams, m);
    end
end
[vLambdaAnaytic, vMultindexToSingleIndexMap] = sort(vLambda_K_before_sort, 'descend');

fprintf(' Before  |  After    |   Multi  | Eigenvalue\n');
fprintf('  sort   |  sort     |   index  | before sort\n');
fprintf('----------------------------------------------\n');
for i = 0:max(sSimParams.PlotSpectM,sSimParams.PlotEigenFuncsM)-1
    m = OneDim2TwoDimIndex(i);
    fprintf('\t%d \t |\t %d\t\t| \t[%d %d]\t|  %f\n', i, vMultindexToSingleIndexMap(i+1)-1,  m(1), m(2), vLambda_K_before_sort(i+1));        
end
end