function res = alg_greedy_mct(P, prm, maxFE) %#ok<INUSD>
%ALG_GREEDY_MCT  List scheduling (LPT order + Minimum Completion Time).
%   Blocks are sorted by decreasing average processing time (LPT rule) and
%   each block is assigned to the miner on which it would finish earliest
%   (MCT rule).  Deterministic, O(n log n + n m).
[n, m] = size(P);
[~, ord] = sort(mean(P, 2), 'descend');
L = zeros(1, m); x = zeros(1, n);
for i = ord'
    [~, j] = min(L + P(i, :));
    x(i) = j; L(j) = L(j) + P(i, j);
end
res.x = x; res.f = max(L); res.fe = 1; res.hist_fe = 1; res.hist_f = res.f;
res.info = struct();
end
