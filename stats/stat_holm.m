function padj = stat_holm(p)
%STAT_HOLM  Holm-Bonferroni adjusted p-values (family-wise error control).
p = p(:)'; k = numel(p);
[ps, o] = sort(p);
adj = zeros(1, k);
for i = 1:k
    adj(i) = min(1, (k - i + 1) * ps(i));
end
for i = 2:k, adj(i) = max(adj(i), adj(i - 1)); end
padj = zeros(1, k); padj(o) = adj;
end
