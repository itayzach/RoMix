function mPhiNumericOrNystrom = FlipSign(sSimParams, mPhiAnalytic, mPhiNumericOrNystrom)
for m = 0:sSimParams.CalcEigenFuncsM-1
    [~, maxAnaIdx] = max(abs(mPhiAnalytic(:,m+1)));
    if sign(mPhiAnalytic(maxAnaIdx,m+1)) ~= sign(mPhiNumericOrNystrom(maxAnaIdx,m+1))
        mPhiNumericOrNystrom(:,m+1) = -mPhiNumericOrNystrom(:,m+1);
    end   
end
end