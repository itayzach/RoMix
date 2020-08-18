function [vLambdaAnalytic, vMultindexToSingleIndexMap, vInd2subRow, vInd2subCol] = CalcAnalyticEigenvalues(sSimParams, sKernelParams, dim, nComponents)

if ~exist('nComponents', 'var')
    nComponents = 1;
end
nEigs = max(sSimParams.PlotSpectM,sSimParams.CalcEigenFuncsM);

if dim == 1
    mLambdaAnaytic = zeros(nComponents, nEigs,1);
    for c = 1:nComponents
        for m = 0:nEigs-1
            mLambdaAnaytic(c,m+1) = lambda(sKernelParams, m);
        end
    end
    
    [ vLambdaAnalytic, vIdx ] = sort(mLambdaAnaytic(:), 'descend');
    [vInd2subRow, vInd2subCol] = ind2sub([nComponents nEigs], vIdx);
    
    vMultindexToSingleIndexMap = 1:nEigs;
elseif dim == 2
    %% Get correct order of eigenvalues (for 1D indexing from multindexing)
    vLambdaBeforeSort = zeros(nEigs,1);
    for i = 0:nEigs-1
        m = OneDim2TwoDimIndex(i);
        vLambdaBeforeSort(i+1) = lambda(sKernelParams, m);
    end
    [vLambdaAnalytic, vMultindexToSingleIndexMap] = sort(vLambdaBeforeSort, 'descend');

    fprintf(' Before  |  After    |   Multi  | Eigenvalue\n');
    fprintf('  sort   |  sort     |   index  | before sort\n');
    fprintf('----------------------------------------------\n');
    for i = 0:nEigs-1
        m = OneDim2TwoDimIndex(i);
        fprintf('\t%d \t |\t %d\t\t| \t[%d %d]\t|  %f\n', i, vMultindexToSingleIndexMap(i+1)-1,  m(1), m(2), vLambdaBeforeSort(i+1));        
    end
    vInd2subRow = [];
    vInd2subCol = [];
end
end