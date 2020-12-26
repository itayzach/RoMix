function mPhiNumericOrNystrom = FlipSign(mPhiAnalytic, mPhiNumericOrNystrom)

M1 = size(mPhiAnalytic,2);
M2 = size(mPhiNumericOrNystrom,2);
M = min(M1,M2);

for m = 0:M-1
    [~, maxAnaIdx] = max(abs(mPhiAnalytic(:,m+1)));
    [~, maxNumIdx] = max(abs(mPhiNumericOrNystrom(:,m+1)));
    if sign(mPhiAnalytic(maxAnaIdx,m+1)) ~= sign(mPhiNumericOrNystrom(maxNumIdx,m+1))
        mPhiNumericOrNystrom(:,m+1) = -mPhiNumericOrNystrom(:,m+1);
    end   
end
end
