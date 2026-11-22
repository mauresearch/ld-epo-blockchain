function draw_box(data, labels, cols)
%DRAW_BOX  Toolbox-free box plot (median, quartiles, 1.5 IQR whiskers,
%   outliers).  data: cell array of vectors.
k = numel(data); hold on;
for i = 1:k
    v = sort(data{i}(:)); q = prctile_local(v, [25 50 75]);
    iqr = q(3) - q(1);
    lo = min(v(v >= q(1) - 1.5 * iqr)); hi = max(v(v <= q(3) + 1.5 * iqr));
    c = cols(i, :);
    patch([i-0.3 i+0.3 i+0.3 i-0.3], [q(1) q(1) q(3) q(3)], c, 'FaceAlpha', 0.35, 'EdgeColor', c);
    plot([i-0.3 i+0.3], [q(2) q(2)], '-', 'Color', c, 'LineWidth', 2);
    plot([i i], [q(3) hi], '-', 'Color', c); plot([i i], [lo q(1)], '-', 'Color', c);
    plot([i-0.12 i+0.12], [hi hi], '-', 'Color', c); plot([i-0.12 i+0.12], [lo lo], '-', 'Color', c);
    out = v(v < lo | v > hi);
    if ~isempty(out), plot(i * ones(size(out)), out, 'o', 'Color', c, 'MarkerSize', 4); end
end
set(gca, 'XTick', 1:k, 'XTickLabel', labels); xlim([0.4 k + 0.6]); box on; hold off;
end

function q = prctile_local(v, p)
n = numel(v); q = zeros(size(p));
for i = 1:numel(p)
    pos = (n - 1) * p(i) / 100 + 1; lo = floor(pos); hi = ceil(pos);
    q(i) = v(lo) + (pos - lo) * (v(hi) - v(lo));
end
end
