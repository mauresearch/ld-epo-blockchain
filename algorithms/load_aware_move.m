function [x, L, moved] = load_aware_move(x, L, P)
%LOAD_AWARE_MOVE  Load-aware neighbourhood operator (the proposed move).
%
%   [x, L, moved] = load_aware_move(x, L, P)
%
%   x : 1 x n assignment vector, L : 1 x m miner loads of x, P : n x m.
%
%   Algorithm:
%     1. Critical (overloaded) miner  j* = argmax_j L(j)   (ties: random)
%     2. TRANSFER neighbourhood: for every block b on j* and every other
%        miner k, the two affected loads become
%            L(j*) - p(b,j*)   and   L(k) + p(b,k)
%        score(b,k) = max of these two values.  The pair with the smallest
%        score is chosen (best-improvement).
%     3. If no transfer gives score < L(j*), the SWAP neighbourhood is
%        examined: exchange b (on j*) with c (on k ~= j*):
%            L(j*) - p(b,j*) + p(c,j*)   and   L(k) - p(c,k) + p(b,k)
%        and the best swap with score < L(j*) is applied.
%     4. If neither neighbourhood improves, x is returned unchanged
%        (moved = false); escaping from this local optimum is left to the
%        exploration phase of D-EPO.
%   The choice is therefore deterministic best-improvement restricted to
%   the critical miner; it has no tunable parameter of its own.  Its
%   application probability (p_LA) and the number of applications per
%   penguin (n_LA) are parameters of LD-EPO.
%   Cost: O(|B_j*| * m) for transfer and O(|B_j*| * n) for swap.

moved = false;
mx = max(L);
cand = find(L == mx);
js = cand(randi(numel(cand)));
B = find(x == js); B = B(:);
if isempty(B), return; end
m = numel(L);

% ---- transfer -------------------------------------------------------
newJs = L(js) - P(B, js);                          % |B| x 1
newK  = bsxfun(@plus, L, P(B, :));                 % |B| x m
sc = max(bsxfun(@times, newJs, ones(1, m)), newK);
sc(:, js) = inf;
[best, idx] = min(sc(:));
if best < L(js)
    [bi, k] = ind2sub(size(sc), idx);
    b = B(bi);
    L(js) = L(js) - P(b, js);
    L(k)  = L(k)  + P(b, k);
    x(b) = k;
    moved = true;
    return;
end

% ---- swap ------------------------------------------------------------
C = find(x ~= js); C = C(:)';
if isempty(C), return; end
kc = x(C);                                         % miners of the C blocks
pc_own = P(sub2ind(size(P), C(:)', kc(:)'));       % 1 x |C|  p(c,k)
pc_js  = P(C, js)';                                % 1 x |C|  p(c,j*)
pb_js  = P(B, js);                                 % |B| x 1  p(b,j*)
pb_k   = P(B, kc);                                 % |B| x |C| p(b,k)
a1 = bsxfun(@plus, L(js) - pb_js, pc_js);          % new load of j*
a2 = bsxfun(@plus, L(kc) - pc_own, pb_k);          % new load of k
sc = max(a1, a2);
[best, idx] = min(sc(:));
if best < L(js)
    [bi, ci] = ind2sub(size(sc), idx);
    b = B(bi); c = C(ci); k = kc(ci);
    L(js) = a1(bi, ci);
    L(k)  = a2(bi, ci);
    x(b) = k; x(c) = js;
    moved = true;
end
end
