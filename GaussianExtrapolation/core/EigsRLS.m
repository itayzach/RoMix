function C = EigsRLS(Phi, gamma1, gamma2, invLambda, L, f, b_maskDataFitTerm)

condNum = cond((Phi).'*(Phi) + gamma1*invLambda);
MTilde = size(Phi,2);
recommendedMTilde = MTilde;
while(condNum > 1e3 && recommendedMTilde >= 10)
    fprintf('Condition number %.3f is too large. Checking with MTilde = %d\n', condNum, recommendedMTilde)
    recommendedMTilde = recommendedMTilde - 10;
    condNum = cond(Phi(1:recommendedMTilde,:));
end
if MTilde ~= recommendedMTilde
    warning('Consider MTilde = %d with condition number %.3f', recommendedMTilde, condNum);
    pause(1);
end

if b_maskDataFitTerm
    assert(size(f,2) == 1, 'The following term would not work for f as a matrix, only as a vector')
    J = diag(sign(abs(f)));
    C = ( (J*Phi).'*(J*Phi) + gamma1*invLambda + gamma2*Phi.'*L*Phi ) \ ((J*Phi).'*f);
else
    C = ( (Phi).'*(Phi) + gamma1*invLambda + gamma2*Phi.'*L*Phi ) \ (Phi.'*f);
end
end