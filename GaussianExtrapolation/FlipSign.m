function [mPhiNumeric, mPhiNystrom] = FlipSign(sSimParams, mPhiAnalytic, mPhiNumeric, mPhiNystrom)
for m = 0:sSimParams.CalcEigenFuncsM-1
    [~, maxAnaIdx] = max(abs(mPhiAnalytic(:,m+1)));
    if sign(mPhiAnalytic(maxAnaIdx,m+1)) ~= sign(mPhiNumeric(maxAnaIdx,m+1))
        mPhiNumeric(:,m+1) = -mPhiNumeric(:,m+1);
    end
    if sign(mPhiAnalytic(maxAnaIdx,m+1)) ~= sign(mPhiNystrom(maxAnaIdx,m+1))
        mPhiNystrom(:,m+1) = -mPhiNystrom(:,m+1);
    end    
end
end