function [vLambdaAnalytic, vComponentIndex, vEigIndex] = CalcAnalyticEigenvalues(sSimParams, sKernelParams, dim, nComponents)

nEigs = max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM);

if dim == 1
    mLambdaAnaytic = zeros(nComponents, nEigs);
    for c = 1:nComponents
        for m = 0:nEigs-1
            mLambdaAnaytic(c,m+1) = lambda(sKernelParams, c, m);
        end
    end
    
    [ vLambdaAnalytic, vIdx ] = sort(mLambdaAnaytic(:), 'descend');
    [vComponentIndex, vEigIndex] = ind2sub([nComponents nEigs], vIdx);
    
elseif dim == 2
    %% Get correct order of eigenvalues (for 1D indexing from multindexing)
    mLambdaBeforeSort = zeros(nComponents, nEigs);
    for c = 1:nComponents
        for i = 0:nEigs-1
            m = OneDim2TwoDimIndex(i);
            mLambdaBeforeSort(c,i+1) = lambda(sKernelParams, c, m);
        end
    end
    vLambdaBeforeSort = mLambdaBeforeSort(:);
    [vLambdaAnalytic, vMultindexToSingleIndexMap] = sort(vLambdaBeforeSort, 'descend');

    fprintf(' Before  |  After   |   Multi   | Eigenvalue  | Eigenvalue  \n');
    fprintf('  sort   |  sort    |   index   | before sort | after sort  \n');
    fprintf('------------------------------------------------------------\n');
    for i = 0:nEigs-1
        m = OneDim2TwoDimIndex(i);
        fprintf('\t%d \t |\t %d\t\t| \t[%d %d]\t|  %f\t  | %f\n', i, vMultindexToSingleIndexMap(i+1)-1,  m(1), m(2), vLambdaBeforeSort(i+1), vLambdaAnalytic(i+1));        
    end
    [vComponentIndex, vEigIndex] = ind2sub([nComponents nEigs], vMultindexToSingleIndexMap);
end
end