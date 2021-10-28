function [mPhi2Flipped, vMinInd] = FlipSign(mPhi1, mPhi2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flip mPhi2 signs and sort them according to mPhi1.
% For each phi1,m:
%     - Align mPhi2 signs according to phi1,m.
%     - Find the k for which the aligned phi2,k is closest to phi1,m 
%       and return mPhi2Flipped(:,m) = mPhi2AlignedToPhi1(:,minInd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[N1, M1] = size(mPhi1);
[N2, M2] = size(mPhi2);
M = min(M1,M2);
assert(N1 == N2, 'lengths must be the same');
vMinInd = zeros(M,1);
mPhi2Flipped = zeros(size(mPhi2));
for m = 1:M
    phi1_m = mPhi1(:,m);
    mPhi2AlignedToPhi1 = nan(size(mPhi2));
    b_vFlippedIndicator = zeros(M,1);
    for k=1:M
        [mPhi2AlignedToPhi1(:,k), b_vFlippedIndicator(k)] = FlipTwoByNormDiff(phi1_m, mPhi2(:,k));
    end
    vNorms = zeros(M,1);
    for k=1:M
        vNorms(k) = norm(phi1_m - mPhi2AlignedToPhi1(:,k));
    end
    [~,minInd] = min(vNorms);
    if minInd ~= m
        warning('m = %d, minInd = %d', m, minInd)
    end
    vMinInd(m) = minInd;
    mPhi2Flipped(:,m) = mPhi2AlignedToPhi1(:,minInd);
end


end

function [phi2Flipped, b_flippedIndicator] = FlipTwoByNormDiff(phi1,phi2)
    normSum = norm(phi1 + phi2);
    normDiff = norm(phi1 - phi2);
    if (normDiff > normSum) && normSum > 1
        phi2Flipped = -phi2;
        b_flippedIndicator = true;
    else
        phi2Flipped = phi2;
        b_flippedIndicator = false;
    end
end
