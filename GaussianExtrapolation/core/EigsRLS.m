function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

lastwarn(''); % clear last warning
if b_maskDataFitTerm
    assert(size(f,2) == 1, 'The following term would not work for f as a matrix, only as a vector')
    J = diag(sign(abs(f)));
    PhiMasked = J*Phi;
    matForInv = PhiMasked.'*PhiMasked + gamma1*invLambda + gamma2*Phi.'*L*Phi;
    C = matForInv \ (PhiMasked.'*f);
else
    matForInv = Phi.'*Phi + gamma1*invLambda + gamma2*Phi.'*L*Phi;
    C = matForInv \ (Phi.'*f);
end
warnMsg = lastwarn();
if ~isempty(warnMsg) % catch warning
    error('\n returned warning: %s\n', warnMsg);
end
fprintf('*********************************************************\n');
fprintf('condition number = %.2f\n', cond(matForInv));
fprintf('*********************************************************\n');
end