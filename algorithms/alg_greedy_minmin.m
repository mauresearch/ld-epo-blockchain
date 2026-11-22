function res = alg_greedy_minmin(P, prm, maxFE) %#ok<INUSD>
%ALG_GREEDY_MINMIN  Min-Min list-scheduling heuristic (Ibarra & Kim 1977;
%   Braun et al. 2001).  At every step, for each unassigned block the
%   minimum completion time over all miners is computed, and the block with
%   the smallest such value is assigned to its best miner.  Deterministic,
%   O(n^2 m).
[n, m] = size(P);
L = zeros(1, m); x = zeros(1, n);
U = true(n, 1);
for step = 1:n
    CT = bsxfun(@plus, P, L);            % completion times n x m
    CT(~U, :) = inf;
    [ctmin, jmin] = min(CT, [], 2);
    [~, i] = min(ctmin);
    j = jmin(i);
    x(i) = j; L(j) = L(j) + P(i, j); U(i) = false;
end
res.x = x; res.f = max(L); res.fe = 1; res.hist_fe = 1; res.hist_f = res.f;
res.info = struct();
end
