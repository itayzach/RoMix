function [vLambdaAnalytic, vComponentIndex, vEigIndex] = CalcAnalyticEigenvalues(nEigs, sKernelParams)

fprintf('Calculating %d eigenvalues... ', nEigs)
nComponents = sKernelParams.sDistParams.estNumComponents;
dim = sKernelParams.sDistParams.dim;

mLambdaBeforeSort = zeros(nComponents, nEigs);
for c = 1:nComponents
    for i = 0:nEigs-1
        m = OneDim2MultiDimIndex(i,dim);
        mLambdaBeforeSort(c,i+1) = lambda(sKernelParams, c, m);
    end
end
vLambdaBeforeSort = mLambdaBeforeSort(:);
[vLambdaAnalytic, vMultindexToSingleIndexMap] = sort(vLambdaBeforeSort, 'descend');

% fprintf('\n')
% fprintf(' Before  |  After   |   Multi   | Eigenvalue  | Eigenvalue  \n');
% fprintf('  sort   |  sort    |   index   | before sort | after sort  \n');
% fprintf('------------------------------------------------------------\n');
% for i = 0:nEigs-1
%     m = OneDim2MultiDimIndex(i,dim);
%     fprintf('\t%3d \t |\t %3d\t\t| \t[%2d %2d]\t|  %8.3f\t  | %8.3f\n', i, vMultindexToSingleIndexMap(i+1)-1,  m(1), m(2), vLambdaBeforeSort(i+1), vLambdaAnalytic(i+1));        
% end

[vComponentIndex, vEigIndex] = ind2sub([nComponents nEigs], vMultindexToSingleIndexMap);
fprintf('Done.\n')
end