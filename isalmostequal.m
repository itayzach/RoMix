function [] = isalmostequal(A, B)
assert(all(all(A - B < 1e-13)));
end