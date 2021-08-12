function [] = isalmostequal(A, B, threshold, msg, do_assert)

if ~exist('threshold', 'var')
    threshold = 1e-16;
end

if ~exist('do_assert', 'var')
    do_assert = true;
end
errormsg = ['You got max error of: ' num2str(max(max(abs(A - B)))) '. Threshold = ' num2str(threshold)];
if do_assert
    if exist('msg', 'var')
        assert(all(all(abs(A - B) <= threshold)), [newline msg newline newline errormsg]);
    else
        assert(all(all(abs(A - B) <= threshold)), [newline newline errormsg]);
    end
else
    fprintf([newline newline errormsg '\n'])
end
end