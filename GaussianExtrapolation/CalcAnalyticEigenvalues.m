function [vLambdaAnaytic, vMultindexToSingleIndexMap] = CalcAnalyticEigenvalues(sSimParams, sKernelParams, dim)

if dim == 1
    vLambdaAnaytic = zeros(max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM),1);
    for m = 0:max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM)-1
        vLambdaAnaytic(m+1) = lambda(sKernelParams, m);
    end
    vMultindexToSingleIndexMap = 1:max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM);
elseif dim == 2
    %% Get correct order of eigenvalues (for 1D indexing from multindexing)
    vLambda_K_before_sort = zeros(max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM),1);
    for i = 0:max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM)-1
        m = OneDim2TwoDimIndex(i);
        vLambda_K_before_sort(i+1) = lambda(sKernelParams, m);
    end
    [vLambdaAnaytic, vMultindexToSingleIndexMap] = sort(vLambda_K_before_sort, 'descend');

    fprintf(' Before  |  After    |   Multi  | Eigenvalue\n');
    fprintf('  sort   |  sort     |   index  | before sort\n');
    fprintf('----------------------------------------------\n');
    for i = 0:max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM)-1
        m = OneDim2TwoDimIndex(i);
        fprintf('\t%d \t |\t %d\t\t| \t[%d %d]\t|  %f\n', i, vMultindexToSingleIndexMap(i+1)-1,  m(1), m(2), vLambda_K_before_sort(i+1));        
    end
end
end