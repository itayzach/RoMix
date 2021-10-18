function [] = isalmostequal(A, B, threshold, msg, do_assert)

if ~exist('threshold', 'var')
    threshold = 1e-16;
end

if ~exist('do_assert', 'var')
    do_assert = true;
end
[maxVal, maxInd] = max(max(abs(A - B)));
[maxRow, maxCol] = ind2sub(size(A), maxInd);
errormsg = sprintf('You got max error of: %d. Threshold = %d.\nA(%d,%d) = %d\nB(%d,%d) = %d\n', ...
   maxVal, threshold, maxRow, maxCol, A(maxRow, maxCol), maxRow, maxCol, B(maxRow, maxCol));
if all(all(abs(A - B) > threshold))
    if do_assert
        if exist('msg', 'var')
            error([newline msg newline newline errormsg]);
        else
            error([newline newline errormsg]);
        end
    else
        fprintf([newline newline errormsg newline newline])
    end
end
end