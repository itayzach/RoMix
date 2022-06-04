function [vLambdaAnalytic, vComponentIndex, vEigIndex, t] = CalcAnalyticEigenvalues(nEigs, sKernelParams)

fprintf('Calculating %d eigenvalues... ', nEigs)
ts = tic;
nComponents = sKernelParams.sDistParams.estNumComponents;
dim = sKernelParams.sDistParams.dim;

OneDim2MultiDimIndexMatrix = LoadOneDim2MultiDimIndexMatrix(nEigs,dim);
mLambdaBeforeSort = zeros(nComponents, nEigs);
for c = 1:nComponents
    p = sKernelParams.sDistParams.GMModel.ComponentProportion(c);
    for i = 1:nEigs
        m = OneDim2MultiDimIndexMatrix(i,:);
        mLambdaBeforeSort(c,i) = p*lambdaD(sKernelParams, c, m);
    end
end
vLambdaBeforeSort = mLambdaBeforeSort(:);
[vLambdaAnalytic, vMultindexToSingleIndexMap] = sort(vLambdaBeforeSort, 'descend');
vLambdaAnalytic = vLambdaAnalytic(1:nEigs);
vMultindexToSingleIndexMap = vMultindexToSingleIndexMap(1:nEigs);

% fprintf('\n')
% fprintf(' Before  |  After   |   Multi   | Eigenvalue  | Eigenvalue  \n');
% fprintf('  sort   |  sort    |   index   | before sort | after sort  \n');
% fprintf('------------------------------------------------------------\n');
% for i = 1:nEigs
%     m = OneDim2MultiDimIndexMatrix(i,:);
%     fprintf('\t%3d \t |\t %3d\t\t|\t[', i-1, vMultindexToSingleIndexMap(i)-1);
%     fprintf('%2d ', m);
%     fprintf(']\t|  %8.3f\t  | %8.3f\n', vLambdaBeforeSort(i), vLambdaAnalytic(i));
% end

[vComponentIndex, vEigIndex] = ind2sub([nComponents nEigs], vMultindexToSingleIndexMap);
t = toc(ts);
fprintf('Done (took %.2f sec).\n', t)
end