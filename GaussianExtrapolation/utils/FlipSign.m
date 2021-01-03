function mPhiNumericOrNystrom = FlipSign(mPhiAnalytic, mPhiNumericOrNystrom)

[N1, M1] = size(mPhiAnalytic);
[N2, M2] = size(mPhiNumericOrNystrom);
M = min(M1,M2);

for m = 0:M-1
    [~, maxAnaIdx] = max(abs(mPhiAnalytic(:,m+1)));
    [~, maxNumIdx] = max(abs(mPhiNumericOrNystrom(:,m+1)));
    if N1 == N2
        if norm(mPhiAnalytic(:,m+1) - mPhiNumericOrNystrom(:,m+1)) >  norm(mPhiAnalytic(:,m+1) + mPhiNumericOrNystrom(:,m+1))
            mPhiNumericOrNystrom(:,m+1) = -mPhiNumericOrNystrom(:,m+1);
        end
    elseif sign(mPhiAnalytic(maxAnaIdx,m+1)) ~= sign(mPhiNumericOrNystrom(maxNumIdx,m+1))
        mPhiNumericOrNystrom(:,m+1) = -mPhiNumericOrNystrom(:,m+1);
    end
end
end
