function [f, L] = eval_pop(X, P)
%EVAL_POP  Makespan of every row of X (vectorised).
%   X : N x n integer matrix, X(r,i) = miner of block i in solution r
%   P : n x m processing-time matrix
%   f : N x 1 makespan  (Fit(X) = max_j sum_i x_ij p_ij)
%   L : N x m completion time of each miner (T_j)
[N, n] = size(X);
m = size(P, 2);
cols = repmat(1:n, N, 1);
pv = P(sub2ind(size(P), cols, X));          % p_{i, X(r,i)}
rows = repmat((1:N)', 1, n);
L = accumarray([rows(:), X(:)], pv(:), [N, m]);
f = max(L, [], 2);
end
