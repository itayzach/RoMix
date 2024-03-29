function [mPhi2Flipped, b_vFlippedIndicator, vMinInd] = FlipSign(mPhi1Base, mPhi2Candidate, b_pairwiseFlip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% b_pairwiseFlipSign = true  --> Navie pairwise flip:
%   For each pair m:
%       - if norm(phi1,m + phi2,m) < norm(phi1,m - phi2,m):
%             return -phi2,m
%         else
%             return phi2,m
% b_pairwiseFlipSign = false --> Flip mPhi2 signs and sort them according to mPhi1.
%   For each phi1,m:
%       - Align mPhi2 signs according to phi1,m.
%       - Find the k for which the aligned phi2,k is closest to phi1,m 
%           and return mPhi2Flipped(:,m) = mPhi2AlignedToPhi1(:,minInd)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('b_pairwiseFlip', 'var')
    b_pairwiseFlip = true;
end
[N1, M1] = size(mPhi1Base);
[N2, M2] = size(mPhi2Candidate);
M = min(M1,M2);
assert(N1 == N2, 'lengths must be the same');
vMinInd = zeros(M,1);
mPhi2Flipped = zeros(size(mPhi2Candidate));

if b_pairwiseFlip
    b_vFlippedIndicator = zeros(M,1);
    for m = 1:M
        [mPhi2Flipped(:,m), b_vFlippedIndicator(m)] = FlipTwoByNormDiff(mPhi1Base(:,m), mPhi2Candidate(:,m));
    end   
else
    for m = 1:M
        phi1_m = mPhi1Base(:,m);
        mPhi2AlignedToPhi1 = nan(size(mPhi2Candidate));
        b_vFlippedIndicator = zeros(M,1);
        for k=1:M
            [mPhi2AlignedToPhi1(:,k), b_vFlippedIndicator(k)] = FlipTwoByNormDiff(phi1_m, mPhi2Candidate(:,k));
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

end

function [phi2Flipped, b_flippedIndicator] = FlipTwoByNormDiff(phi1,phi2)
    normSum = norm(phi1 + phi2);
    normDiff = norm(phi1 - phi2);
    if (normDiff > normSum)
        phi2Flipped = -phi2;
        b_flippedIndicator = true;
    else
        phi2Flipped = phi2;
        b_flippedIndicator = false;
    end
end
