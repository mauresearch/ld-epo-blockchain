function draw_gantt(x, P, ttl, show_labels)
%DRAW_GANTT  Gantt chart of an assignment: one row per miner (ALL m miners
%   are drawn, including idle ones), blocks in index order.
if nargin < 4, show_labels = true; end
[n, m] = size(P);
cols = lines(7); cols = [cols; 0.6 0.8 0.6; 0.9 0.7 0.4; 0.7 0.7 0.9];
hold on;
for j = 1:m
    B = find(x == j); t = 0;
    for b = B
        d = P(b, j);
        c = cols(mod(b - 1, size(cols, 1)) + 1, :);
        patch([t t+d t+d t], [j-0.4 j-0.4 j+0.4 j+0.4], c, 'EdgeColor', 'k', 'FaceAlpha', 0.85);
        if show_labels
            text(t + d / 2, j, sprintf('B%d', b), 'HorizontalAlignment', 'center', ...
                 'FontSize', 7);
        end
        t = t + d;
    end
end
L = accumarray(x(:), P(sub2ind(size(P), (1:n)', x(:))), [m 1]);
Cmax = max(L);
plot([Cmax Cmax], [0.4 m + 0.6], 'r--', 'LineWidth', 1.5);
text(Cmax, m + 0.75, sprintf('C_{max} = %g', Cmax), 'Color', 'r', 'HorizontalAlignment', 'right');
set(gca, 'YTick', 1:m, 'YTickLabel', arrayfun(@(j) sprintf('M%d', j), 1:m, 'UniformOutput', false));
ylim([0.4 m + 1]); xlim([0 Cmax * 1.05]);
xlabel('Time'); ylabel('Miner'); title(ttl); box on; hold off;
end
