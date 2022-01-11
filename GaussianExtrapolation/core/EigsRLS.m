function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

lastwarn(''); % clear last warning
if b_maskDataFitTerm
    if(size(f,2) == 10)
        assert(islogical(f), 'If you got here, your f must be one hot encodeing classes matrix')
        f_sum = sum(f,2);
        J = diag(f_sum);
    elseif(size(f,2) == 1)
        J = diag(abs(f));
    else
        error('Semi-Supervised Learning masking is avaliable only for USPS and TwoMoons')
    end
    JPhi = J*Phi;
else
    JPhi = Phi;
end
if gamma2 == 0 || isempty(L)
    matToInv = JPhi.'*JPhi + gamma1*invLambda;
else
    matToInv = JPhi.'*JPhi + gamma1*invLambda + gamma2*Phi.'*L*Phi;
end
C = matToInv \ (JPhi.'*f);

warnMsg = lastwarn();
if ~isempty(warnMsg) % catch warning
    error('\n returned warning: %s\n', warnMsg);
end
fprintf('*********************************************************\n');
fprintf('cond(matToInv) = %.2f\n', cond(matToInv));
fprintf('*********************************************************\n');
end