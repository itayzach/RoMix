function [vLabeledInd, J] = GetUnlabeledNodesMask(f)
if(size(f,2) == 10) % 10 class dataset - MNIST/USPS
    assert(isequal(f,sign(f)), 'If you got here, your f must be one hot encodeing classes matrix')
    f_sum = sum(f,2);
    J = diag(f_sum);
    vLabeledInd = f_sum;
else
    % Build J according to first signal
    vLabeledInd = abs(f(:,1)) > 0;
    J = diag(vLabeledInd);

    % Now make sure all signals are labeled/unlabled in the same entries before we continue
    mIsLabeled = abs(f) > 0;
    mIsFirstColDup = repmat(abs(f(:,1)) > 0,1,size(f,2));
    assert(isequal(mIsFirstColDup,mIsLabeled))
end
vLabeledInd = find(vLabeledInd);
end