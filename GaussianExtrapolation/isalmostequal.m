function [] = isalmostequal(A, B, threshold, msg)
if exist('msg', 'var')
    assert(all(all(A - B < threshold)), [newline msg newline newline 'You got error of: ' num2str(abs(sum(sum(A - B))))]);
else
    assert(all(all(A - B < threshold)));
end
end