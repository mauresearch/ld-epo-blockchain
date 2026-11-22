function [p, chi2, avgrank, pID, FF] = stat_friedman(X)
%STAT_FRIEDMAN  Friedman test over a (blocks x methods) matrix X
%   (rows = problem instances, columns = methods, smaller is better).
%   Also returns the Iman-Davenport correction (FF, pID).
%   Demsar (2006), Derrac et al. (2011).
[b, k] = size(X);
R = zeros(b, k);
for i = 1:b
    R(i, :) = rank_row(X(i, :));
end
avgrank = mean(R, 1);
chi2 = 12 * b / (k * (k + 1)) * (sum(avgrank.^2) - k * (k + 1)^2 / 4);
p = 1 - gammainc(chi2 / 2, (k - 1) / 2);
FF = (b - 1) * chi2 / (b * (k - 1) - chi2);
d1 = k - 1; d2 = (k - 1) * (b - 1);
if FF > 0 && isfinite(FF)
    pID = 1 - betainc(d1 * FF / (d1 * FF + d2), d1 / 2, d2 / 2);
else
    pID = 1;
end
end

function r = rank_row(x)
[xs, ord] = sort(x); n = numel(x); r = zeros(1, n); i = 1;
while i <= n
    j = i;
    while j < n && xs(j + 1) == xs(i), j = j + 1; end
    r(ord(i:j)) = (i + j) / 2; i = j + 1;
end
end
