function [C, t] = RoMix(Phi, gamma1, gamma2, lambdaPhi, L, f)

ts = tic;
invLambda = diag(1./lambdaPhi);

vLabeledInd = GetUnlabeledNodesMask(f);
JPhi = Phi(vLabeledInd,:);

lastwarn(''); % clear last warning
if gamma2 == 0 || isempty(L)
    matToInv = JPhi.'*JPhi + gamma1*invLambda;
else
    assert(gamma2 > 0 && ~isempty(L))
    matToInv = JPhi.'*JPhi + gamma1*invLambda + gamma2*Phi.'*L*Phi;
end
C = matToInv \ (JPhi.'*f(vLabeledInd,:));

warnMsg = lastwarn();
if ~isempty(warnMsg) % catch warning
    error('\n returned warning: %s\n', warnMsg);
end
t = toc(ts);
PrintMethodStats('RoMix',Phi,C,invLambda,L);

end
