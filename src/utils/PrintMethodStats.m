function PrintMethodStats(methodName,Phi,C,invLambda,L,matToInv)
fRec = Phi*C;
nFuncs = size(fRec,2);

if ~isempty(invLambda)
    fprintf('%s: c^T*invLambda*c = \n\t', methodName)
    for i = 1:nFuncs
        fprintf('\t(%2d)%8.2f ',i,C(:,i).'*invLambda*C(:,i))
        if mod(i,10) == 0, fprintf('\n\t'), end
    end
    fprintf('\n')
end

if ~isempty(L)
    fprintf('%s: f^T*L*f = \n\t', methodName)
    for i = 1:nFuncs
        fprintf('\t(%2d)%8.2f ',i,fRec(:,i).'*L*fRec(:,i))
        if mod(i,10) == 0, fprintf('\n\t'), end
    end
    fprintf('\n')
end

fprintf('*********************************************************\n');
if exist('matToInv', 'var') && ~isempty(matToInv)
    fprintf('* cond(matToInv) = %.2f\n', cond(matToInv));
    fprintf('*********************************************************\n');
end
end
