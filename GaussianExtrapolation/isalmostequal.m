function [] = isalmostequal(A, B, threshold, msg)
if exist('msg', 'var')
    assert(all(all(A - B < threshold)), msg);
else
    assert(all(all(A - B < threshold)));
end
end