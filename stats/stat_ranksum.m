function [p, z] = stat_ranksum(a, b)
%STAT_RANKSUM  Two-sided Wilcoxon rank-sum (Mann-Whitney U) test with tie
%   correction and normal approximation (no toolbox needed; equivalent to
%   ranksum(a,b,'method','approximate') of the Statistics Toolbox).
a = a(:); b = b(:);
n1 = numel(a); n2 = numel(b); N = n1 + n2;
[r, t] = tiedrank_local([a; b]);
W = sum(r(1:n1));
mu = n1 * (N + 1) / 2;
sig2 = n1 * n2 / 12 * ((N + 1) - t / (N * (N - 1)));
if sig2 <= 0, p = 1; z = 0; return; end
z = (W - mu - 0.5 * sign(W - mu)) / sqrt(sig2);       % continuity correction
p = erfc(abs(z) / sqrt(2));
p = min(1, p);
end

function [r, tieadj] = tiedrank_local(x)
[xs, ord] = sort(x);
n = numel(x); r = zeros(n, 1); tieadj = 0;
i = 1;
while i <= n
    j = i;
    while j < n && xs(j + 1) == xs(i), j = j + 1; end
    r(ord(i:j)) = (i + j) / 2;
    tcount = j - i + 1;
    tieadj = tieadj + (tcount^3 - tcount);
    i = j + 1;
end
end
