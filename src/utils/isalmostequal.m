function [] = isalmostequal(A, B, threshold, msg, b_assertOnFail, b_printMaxErr)

if ~exist('threshold', 'var')
    threshold = 1e-16;
end

if ~exist('b_printMaxErr', 'var')
    b_printMaxErr = false;
end

if ~exist('b_assertOnFail', 'var')
    b_assertOnFail = true;
end
[maxVal, maxInd] = max(max(abs(A - B)));
[maxRow, maxCol] = ind2sub(size(A), maxInd);
errormsg = sprintf('Max error = %d \nThreshold = %d\nA(%d,%d) = %d\nB(%d,%d) = %d\n', ...
   maxVal, threshold, maxRow, maxCol, A(maxRow, maxCol), maxRow, maxCol, B(maxRow, maxCol));
if any(abs(A(:) - B(:)) > threshold)
    if b_assertOnFail
        if exist('msg', 'var')
            error([newline msg newline newline errormsg]);
        else
            error([newline newline errormsg]);
        end
    else
        fprintf([newline newline errormsg newline newline])
    end
end
if b_printMaxErr
    fprintf(errormsg);
end
end