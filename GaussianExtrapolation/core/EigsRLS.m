function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

if b_maskDataFitTerm
    J = GetUnlabeledNodesMask(f);
else
    J = eye(size(f,1));
end
JPhi = J*Phi;

lastwarn(''); % clear last warning
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

PrintEigsRLSStats(Phi,C,invLambda,L,matToInv)

end
