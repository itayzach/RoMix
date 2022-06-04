function [alpha, t] = LapRLS(W,f,L,lambda1,lambda2,interpRatio,b_renormalize, b_maskDataFitTerm)

n=size(W,1); % total examples
I=eye(n);

if b_maskDataFitTerm
    J = GetUnlabeledNodesMask(f);
else
    J = I;
end

lastwarn(''); % clear last warning
ts = tic;
if ~isempty(L) && lambda2~=0
    matToInv = J*W + lambda1*I + lambda2*L*W;
else
    matToInv = J*W + lambda1*I;
end
alpha = matToInv\f;
if b_renormalize
    alpha = alpha/sqrt(interpRatio);
end
t = toc(ts);

PrintMethodStats('LapRLS',W,alpha,[],L);
end