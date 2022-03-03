function PrintEigsRLSStats(Phi,C,invLambda,L,matToInv)
fRec = Phi*C;
nFuns = size(fRec,2);
fprintf('EigsRLS: c^T*invLambda*c = \n\t')
for i = 1:nFuns
    fprintf('\t(%2d)%8.2f ',i,C(:,i).'*invLambda*C(:,i))
    if mod(i,10) == 0, fprintf('\n\t'), end
end
fprintf('\n')

if ~isempty(L)
    fprintf('EigsRLS: f^T*L*f = \n\t')
    for i = 1:nFuns
        fprintf('\t(%2d)%8.2f ',i,fRec(:,i).'*L*fRec(:,i))
        if mod(i,10) == 0, fprintf('\n\t'), end
    end
    fprintf('\n')
end

fprintf('*********************************************************\n');
fprintf('* cond(matToInv) = %.2f\n', cond(matToInv));
fprintf('*********************************************************\n');
end
