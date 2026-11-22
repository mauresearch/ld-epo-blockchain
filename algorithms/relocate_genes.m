function [y, Ly] = relocate_genes(y, Ly, P, k)
%RELOCATE_GENES  Exploration (relocation) operator of LD-EPO.
%   k randomly chosen blocks are moved to a randomly chosen miner among the
%   miners that can receive them WITHOUT exceeding the current makespan
%   (plateau / non-worsening relocation).  If no such miner exists the
%   target miner is chosen uniformly at random.  The operator diversifies
%   the load profile of the non-critical miners, which in turn opens new
%   improving transfers/swaps for the load-aware move.
[n, m] = size(P);
Cmax = max(Ly);
pos = randperm(n, min(k, n));
for b = pos
    cur = y(b);
    newL = Ly + P(b, :);
    ok = find(newL <= Cmax);
    ok(ok == cur) = [];
    if isempty(ok)
        j = randi(m);
    else
        j = ok(randi(numel(ok)));
    end
    if j ~= cur
        Ly(cur) = Ly(cur) - P(b, cur);
        Ly(j) = Ly(j) + P(b, j);
        y(b) = j;
    end
end
end
