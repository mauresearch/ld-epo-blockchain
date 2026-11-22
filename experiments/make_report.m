function make_report()
%MAKE_REPORT  Builds every table (CSV) and figure (PNG) of the paper from the
%   saved raw results.
cfg = config_experiment();
od = cfg.outdir; fd = cfg.figdir;
if ~exist(fd, 'dir'), mkdir(fd); end
set(0, 'defaultfigurevisible', 'off');

% ---------------- load ---------------------------------------------------
RS = {};
for s = 1:numel(cfg.scen)
    fr = fullfile(od, sprintf('raw_%s.mat', cfg.scen(s).name));
    if exist(fr, 'file'), S = load(fr); RS{end+1} = S.res; end %#ok<AGROW>
end
if isempty(RS), error('no raw results found - run run_experiments first'); end
nS = numel(RS);
mnames = {RS{1}.methods.name};
nM = numel(mnames);
iP = find(strcmp(mnames, 'LD-EPO'));

% MILP reference (optional)
opt = nan(1, nS); optlb = nan(1, nS); optst = repmat({''}, 1, nS);
fm = fullfile(od, 'milp_results.csv');
if exist(fm, 'file')
    fid = fopen(fm); fgetl(fid);
    while true
        ln = fgetl(fid); if ~ischar(ln), break; end
        c = strsplit(ln, ',');
        k = find(cellfun(@(r) strcmp(r.scen.name, c{1}), RS));
        if ~isempty(k)
            opt(k) = str2double(c{4}); optlb(k) = str2double(c{5}); optst{k} = c{6};
        end
    end
    fclose(fid);
end

% ---------------- Table: main results -----------------------------------
fid = fopen(fullfile(od, 'table_main.csv'), 'w');
fprintf(fid, 'scenario,m,n,method,role,min,max,mean,std,median,time_s,FE,RPD_best_known_%%,gap_to_MILP_LB_%%\n');
MEAN = zeros(nS, nM);
for s = 1:nS
    R = RS{s}; M = R.methods;
    bk = min(cellfun(@min, {M.f}));                 % best-known value
    if ~isnan(opt(s)), bk = min(bk, opt(s)); end
    for k = 1:nM
        f = M(k).f; MEAN(s, k) = mean(f);
        gap = nan; if ~isnan(optlb(s)), gap = 100 * (mean(f) - optlb(s)) / optlb(s); end
        fprintf(fid, '%s,%d,%d,%s,%s,%g,%g,%.3f,%.3f,%g,%.3f,%d,%.3f,%.3f\n', R.scen.name, ...
            R.scen.m, R.scen.n, M(k).name, M(k).role, min(f), max(f), mean(f), std(f), ...
            median(f), mean(M(k).time), round(mean(M(k).fe)), 100 * (mean(f) - bk) / bk, gap);
    end
end
fclose(fid);

% ---------------- Table: statistical tests vs LD-EPO --------------------
fid = fopen(fullfile(od, 'table_stats.csv'), 'w');
fprintf(fid, 'scenario,method,mean_method,mean_LD-EPO,p_wilcoxon,p_holm,A12_LD-EPO_better,result\n');
WTL = zeros(nM, 3);                                  % win / tie / loss of LD-EPO
for s = 1:nS
    M = RS{s}.methods; fp = M(iP).f;
    others = setdiff(1:nM, iP);
    p = ones(1, numel(others)); A = p;
    for q = 1:numel(others)
        p(q) = stat_ranksum(fp, M(others(q)).f);
        A(q) = stat_a12(fp, M(others(q)).f);
    end
    ph = stat_holm(p);
    for q = 1:numel(others)
        k = others(q);
        if ph(q) < cfg.alpha && mean(fp) < mean(M(k).f), r = '+'; WTL(k, 1) = WTL(k, 1) + 1;
        elseif ph(q) < cfg.alpha && mean(fp) > mean(M(k).f), r = '-'; WTL(k, 3) = WTL(k, 3) + 1;
        else, r = '='; WTL(k, 2) = WTL(k, 2) + 1; end
        fprintf(fid, '%s,%s,%.3f,%.3f,%.3e,%.3e,%.3f,%s\n', RS{s}.scen.name, M(k).name, ...
            mean(M(k).f), mean(fp), p(q), ph(q), A(q), r);
    end
end
fclose(fid);
fid = fopen(fullfile(od, 'table_win_tie_loss.csv'), 'w');
fprintf(fid, 'method,LD-EPO_better(+),no_significant_difference(=),LD-EPO_worse(-)\n');
for k = setdiff(1:nM, iP)
    fprintf(fid, '%s,%d,%d,%d\n', mnames{k}, WTL(k, 1), WTL(k, 2), WTL(k, 3));
end
fclose(fid);

% ---------------- Table: Friedman ----------------------------------------
fid = fopen(fullfile(od, 'table_friedman.csv'), 'w');
[pF, chi2, ar, pID, FF] = stat_friedman(MEAN);
[~, o] = sort(ar);
fprintf(fid, 'rank_position,method,average_rank\n');
for i = 1:nM, fprintf(fid, '%d,%s,%.3f\n', i, mnames{o(i)}, ar(o(i))); end
fprintf(fid, '\nFriedman chi2,%.4f\np-value,%.4e\nIman-Davenport F,%.4f\np-value (ID),%.4e\nblocks (scenarios),%d\n', ...
    chi2, pF, FF, pID, nS);
fclose(fid);

% ---------------- Table: improvement (%) of LD-EPO -----------------------
fid = fopen(fullfile(od, 'table_improvement.csv'), 'w');
fprintf(fid, 'scenario');
for k = setdiff(1:nM, iP), fprintf(fid, ',vs_%s', mnames{k}); end
fprintf(fid, '\n');
for s = 1:nS
    fprintf(fid, '%s', RS{s}.scen.name);
    for k = setdiff(1:nM, iP)
        fprintf(fid, ',%.2f', 100 * (MEAN(s, k) - MEAN(s, iP)) / MEAN(s, k));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% ---------------- Table: MILP gap ----------------------------------------
if any(~isnan(opt))
    fid = fopen(fullfile(od, 'table_milp_gap.csv'), 'w');
    fprintf(fid, 'scenario,MILP_best,MILP_lower_bound,status,LD-EPO_best,LD-EPO_mean,gap_best_%%,gap_mean_%%,runs_reaching_MILP_best\n');
    for s = find(~isnan(opt))
        f = RS{s}.methods(iP).f;
        fprintf(fid, '%s,%g,%g,%s,%g,%.3f,%.3f,%.3f,%d/%d\n', RS{s}.scen.name, opt(s), optlb(s), ...
            optst{s}, min(f), mean(f), 100 * (min(f) - optlb(s)) / optlb(s), ...
            100 * (mean(f) - optlb(s)) / optlb(s), sum(f <= opt(s)), numel(f));
    end
    fclose(fid);
end

% ---------------- Figures ------------------------------------------------
cmap = [0.55 0.55 0.55; 0.20 0.20 0.20; 0.45 0.30 0.15; 0.85 0.33 0.10; 0.93 0.69 0.13; ...
        0.49 0.18 0.56; 0.30 0.75 0.93; 0.64 0.08 0.18; 0.47 0.67 0.19; 0.00 0.45 0.74; ...
        0.75 0.55 0.85; 0.00 0.60 0.30];
stoch = find(~[RS{1}.methods.det] & ~strcmp(mnames, 'Random'));

for s = 1:nS
    R = RS{s}; M = R.methods; nm = R.scen.name;
    % convergence
    h = figure('visible', 'off', 'Position', [100 100 760 480]); hold on;
    hl = []; lg = {};
    for k = stoch
        lw = 1.2; if k == iP, lw = 2.5; end
        hl(end+1) = plot(R.grid, M(k).curve, '-', 'Color', cmap(k, :), 'LineWidth', lw); %#ok<AGROW>
        lg{end+1} = M(k).name; %#ok<AGROW>
    end
    kM = find(strcmp(mnames, 'LPT-MCT'));
    hl(end+1) = plot(R.grid([1 end]), mean(M(kM).f) * [1 1], '--', 'Color', cmap(kM, :), 'LineWidth', 1.2);
    lg{end+1} = 'LPT-MCT';
    if ~isnan(opt(s))
        hl(end+1) = plot(R.grid([1 end]), optlb(s) * [1 1], 'k:', 'LineWidth', 1.5);
        lg{end+1} = 'MILP optimum / lower bound';
    end
    lo = min(cellfun(@min, {M.f})); if ~isnan(optlb(s)), lo = min(lo, optlb(s)); end
    hi = mean(M(kM).f) * 1.12;
    ylim([lo * 0.98, hi]); xlim(R.grid([1 end]));
    xlabel('Function evaluations'); ylabel(sprintf('Mean best makespan (%d runs)', cfg.runs));
    title(sprintf('%s (m = %d, n = %d): convergence', nm, R.scen.m, R.scen.n));
    legend(hl, lg, 'Location', 'northeast', 'FontSize', 8); grid on; box on;
    print(h, '-dpng', '-r150', fullfile(fd, sprintf('fig_convergence_%s.png', nm))); close(h);

    % box plot
    h = figure('visible', 'off', 'Position', [100 100 860 420]);
    draw_box({M(stoch).f}, mnames(stoch), cmap(stoch, :));
    ylabel('Makespan'); title(sprintf('%s: distribution over %d runs', nm, numel(M(iP).f)));
    grid on;
    print(h, '-dpng', '-r150', fullfile(fd, sprintf('fig_boxplot_%s.png', nm))); close(h);

    % Gantt (all miners shown)
    h = figure('visible', 'off', 'Position', [100 100 900 max(300, 28 * R.scen.m + 120)]);
    draw_gantt(M(iP).bestx, R.P, sprintf('%s: best LD-EPO assignment (m = %d miners)', nm, R.scen.m), R.scen.n <= 100);
    print(h, '-dpng', '-r150', fullfile(fd, sprintf('fig_gantt_%s.png', nm))); close(h);
end

% ablation figure: RPD (%) of the mean w.r.t. best-known value
abl = [find(strcmp(mnames, 'EPO-basic')), find(strcmp(mnames, 'D-EPO')), ...
       find(strcmp(mnames, 'LA-LS')), find(strcmp(mnames, 'LD-EPO-RI')), iP];
RPD = zeros(nS, numel(abl));
for s = 1:nS
    bk = min(cellfun(@min, {RS{s}.methods.f})); if ~isnan(opt(s)), bk = min(bk, opt(s)); end
    RPD(s, :) = 100 * (MEAN(s, abl) - bk) / bk;
end
h = figure('visible', 'off', 'Position', [100 100 860 420]);
b = bar(RPD); for i = 1:numel(abl), set(b(i), 'FaceColor', cmap(abl(i), :)); end
set(gca, 'XTickLabel', cellfun(@(r) r.scen.name, RS, 'UniformOutput', false));
ylabel('RPD of mean makespan to best known (%)'); legend(mnames(abl), 'Location', 'northwest');
title('Ablation study'); grid on;
print(h, '-dpng', '-r150', fullfile(fd, 'fig_ablation.png')); close(h);

% runtime figure
h = figure('visible', 'off', 'Position', [100 100 760 460]); hold on;
nn = cellfun(@(r) r.scen.n, RS);
for k = [stoch, find([RS{1}.methods.det])]
    t = cellfun(@(r) mean(r.methods(k).time), RS);
    lw = 1.2; if k == iP, lw = 2.5; end
    plot(nn, t, '-o', 'Color', cmap(k, :), 'LineWidth', lw, 'MarkerSize', 4);
end
set(gca, 'XScale', 'log', 'YScale', 'log');
legend(mnames([stoch, find([RS{1}.methods.det])]), 'Location', 'eastoutside', 'FontSize', 8);
xlabel('Number of blocks n'); ylabel('Mean CPU time per run (s)');
title(sprintf('Running time (budget = %d evaluations)', cfg.maxFE)); grid on; box on;
print(h, '-dpng', '-r150', fullfile(fd, 'fig_runtime.png')); close(h);

% sensitivity figure
fl = dir(fullfile(od, 'sensitivity_*.mat'));
if ~isempty(fl)
    S = struct('param', {}, 'value', {}, 'f', {}, 'time', {});
    for q = 1:numel(fl)
        T = load(fullfile(od, fl(q).name)); S = [S, T.S]; %#ok<AGROW>
    end
    pn = fieldnames(cfg.sens.grid); pn = pn(ismember(pn, unique({S.param})));
    h = figure('visible', 'off', 'Position', [100 100 900 600]);
    for a = 1:numel(pn)
        idx = find(strcmp({S.param}, pn{a}));
        v = [S(idx).value]; mu = arrayfun(@(z) mean(z.f), S(idx)); sd = arrayfun(@(z) std(z.f), S(idx));
        subplot(2, 2, a);
        he = errorbar(1:numel(v), mu, sd, '-o');
        set(he, 'Color', cmap(iP, :), 'LineWidth', 1.5);
        set(gca, 'XTick', 1:numel(v), 'XTickLabel', arrayfun(@(z) num2str(z), v, 'UniformOutput', false));
        xlim([0.5 numel(v) + 0.5]);
        xlabel(strrep(pn{a}, '_', '\_')); ylabel('Mean makespan \pm std'); grid on;
    end
    print(h, '-dpng', '-r150', fullfile(fd, 'fig_sensitivity.png')); close(h);
end
fprintf('Tables written to %s, figures written to %s\n', od, fd);
end
