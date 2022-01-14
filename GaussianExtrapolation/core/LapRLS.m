function alpha = LapRLS(W,f,L,lambda1,lambda2,interpRatio,b_renormalize, b_maskDataFitTerm)

n=size(W,1); % total examples
I=eye(n);

lastwarn(''); % clear last warning
if b_maskDataFitTerm
    if(size(f,2) == 10)
        assert(isequal(f,sign(f)), 'If you got here, your f must be one hot encodeing classes matrix')
        f_sum = sum(f,2);
        J = diag(f_sum);
    elseif(size(f,2) == 1)
        J = diag(abs(f));
    else
        error('Semi-Supervised Learning masking is avaliable only for USPS and TwoMoons')
    end
else
    J = I;
end

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