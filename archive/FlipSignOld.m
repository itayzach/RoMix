function mPhiNumericOrNystrom = FlipSignOld(mPhiAnalytic, mPhiNumericOrNystrom)

[N1, M1] = size(mPhiAnalytic);
[N2, M2] = size(mPhiNumericOrNystrom);
M = min(M1,M2);
assert(N1 == N2, 'lengths must be the same');
for m = 0:M-1
    if sign(mPhiAnalytic(1,m+1)) ~= sign(mPhiNumericOrNystrom(1,m+1))
        mPhiNumericOrNystrom(:,m+1) = -mPhiNumericOrNystrom(:,m+1);
    end
end

end
