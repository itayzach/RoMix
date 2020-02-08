function [] = isalmostequal(A, B, threshold)
assert(all(all(A - B < threshold)));
end