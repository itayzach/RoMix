function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

if b_maskDataFitTerm
    assert(size(f,2) == 1, 'The following term would not work for f as a matrix, only as a vector')
    J = diag(sign(abs(f)));
    C = ( (J*Phi).'*(J*Phi) + gamma1*invLambda + gamma2*Phi.'*L*Phi ) \ ((J*Phi).'*f);
else
    C = ( (Phi).'*(Phi) + gamma1*invLambda + gamma2*Phi.'*L*Phi ) \ (Phi.'*f);
end
end