function [] = isalmostequal(A, B, threshold, msg, do_assert)

if ~exist('do_assert', 'var')
    do_assert = true;
end
if do_assert
    if exist('msg', 'var')
        assert(all(all(A - B <= threshold)), [newline msg newline newline 'You got error of: ' num2str(abs(sum(sum(A - B))))]);
    else
        assert(all(all(A - B <= threshold)), [newline newline 'You got error of: ' num2str(abs(sum(sum(A - B))))]);
    end
else
    fprintf([newline newline 'You got error of: ' num2str(abs(sum(sum(A - B)))) '\n'])
end
end