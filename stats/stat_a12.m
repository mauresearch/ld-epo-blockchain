function A = stat_a12(a, b)
%STAT_A12  Vargha-Delaney A12 effect size for MINIMISATION:
%   probability that a run of method a gives a SMALLER makespan than a run
%   of method b (ties count 1/2).  0.5 = no effect; >0.71 large,
%   >0.64 medium, >0.56 small effect in favour of a.
a = a(:); b = b(:)';
A = (sum(sum(bsxfun(@lt, a, b))) + 0.5 * sum(sum(bsxfun(@eq, a, b)))) / (numel(a) * numel(b));
end
