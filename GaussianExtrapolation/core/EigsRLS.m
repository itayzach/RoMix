function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

lastwarn(''); % clear last warning
if b_maskDataFitTerm
    assert(size(f,2) == 1, 'The following term would not work for f as a matrix, only as a vector')
    J = diag(sign(abs(f)));
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