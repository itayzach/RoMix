function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

lastwarn(''); % clear last warning
if b_maskDataFitTerm
    assert(size(f,2) == 1, 'The following term would not work for f as a matrix, only as a vector')
    J = diag(sign(abs(f)));
    Phi = J*Phi;
end
if gamma2 == 0 || isempty(L)
    matToInv = Phi.'*Phi + gamma1*invLambda;
else
    matToInv = Phi.'*Phi + gamma1*invLambda + gamma2*Phi.'*L*Phi;
end
C = matToInv \ (Phi.'*f);

warnMsg = lastwarn();
if ~isempty(warnMsg) % catch warning
    error('\n returned warning: %s\n', warnMsg);
end
fprintf('*********************************************************\n');
fprintf('cond(matToInv) = %.2f\n', cond(matToInv));
fprintf('*********************************************************\n');
end