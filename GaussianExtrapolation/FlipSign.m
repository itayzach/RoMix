function [mPhiNumeric, mPhiNystrom] = FlipSign(sSimParams, mPhiAnalytic, mPhiNumeric, mPhiNystrom)
for m = 0:sSimParams.CalcEigenFuncsM-1
    % Just in case first entry is zero, find the next non-zero entry
    eigEntry = 1;
    while sign(mPhiAnalytic(eigEntry,m+1)) == 0 || sign(mPhiNumeric(eigEntry,m+1)) == 0
        eigEntry = eigEntry + 1;
    end
    if sign(mPhiAnalytic(eigEntry,m+1)) ~= sign(mPhiNumeric(eigEntry,m+1))
        mPhiNumeric(:,m+1) = -mPhiNumeric(:,m+1);
    end
    
    eigEntry = 1;
    while sign(mPhiAnalytic(eigEntry,m+1)) == 0 || sign(mPhiNystrom(eigEntry,m+1)) == 0
        eigEntry = eigEntry + 1;
    end
    if sign(mPhiAnalytic(eigEntry,m+1)) ~= sign(mPhiNystrom(eigEntry,m+1))
        mPhiNystrom(:,m+1) = -mPhiNystrom(:,m+1);
    end
end
end