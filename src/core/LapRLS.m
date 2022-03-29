function alpha = LapRLS(W,f,L,lambda1,lambda2,interpRatio,b_renormalize, b_maskDataFitTerm)

n=size(W,1); % total examples
I=eye(n);

if b_maskDataFitTerm
    J = GetUnlabeledNodesMask(f);
else
    J = I;
end

lastwarn(''); % clear last warning
if ~isempty(L) && lambda2~=0
    alpha = (J*W + lambda1*I + lambda2*L*W)\f;
else
    alpha = (J*W + lambda1*I)\f;
end
if b_renormalize
    alpha = alpha/sqrt(interpRatio);
end

fRec = W*alpha;
fprintf('LapRLS : first fRec^T*L*fRec = %.4f\n',fRec(:,1).'*L*fRec(:,1))

end