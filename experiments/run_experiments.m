function run_experiments(scen_ids)
%RUN_EXPERIMENTS  Main experiment: every method x every scenario x cfg.runs.
%
%   run_experiments          % all scenarios
%   run_experiments([1 2 3]) % only some (useful to run in parallel)
%
%   Saves results/raw_S<k>.mat containing, for every method, the makespan,
%   CPU time and number of evaluations of each run (needed for the
%   statistical tests), the averaged convergence curve and
%   the best assignment found.  Run seeds are fixed (seed = run index), so
%   every number in the paper can be reproduced exactly.
cfg = config_experiment();
if nargin < 1, scen_ids = 1:numel(cfg.scen); end
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end
grid = linspace(0, cfg.maxFE, cfg.grid);

for s = scen_ids
    sc = cfg.scen(s);
    P = get_instance(cfg, s);
    dlmwrite(fullfile(cfg.outdir, sprintf('P_%s.csv', sc.name)), P);
    fprintf('\n=== %s : m = %d, n = %d ===\n', sc.name, sc.m, sc.n);
    M = struct([]);
    for k = 1:numel(cfg.methods)
        md = cfg.methods(k);
        R = cfg.runs; if md.det, R = 1; end
        f = zeros(1, R); tm = zeros(1, R); fe = zeros(1, R);
        curves = zeros(R, numel(grid));
        bestf = inf; bestx = [];
        for r = 1:R
            rng(r, 'twister');
            t0 = tic;
            o = feval(md.fun, P, md.prm, cfg.maxFE);
            tm(r) = toc(t0);
            fchk = eval_pop(o.x, P);
            if abs(fchk - o.f) > 1e-9
                error('%s returned an inconsistent objective', md.name);
            end
            f(r) = o.f; fe(r) = o.fe;
            curves(r, :) = resample_curve(o.hist_fe, o.hist_f, grid);
            if o.f < bestf, bestf = o.f; bestx = o.x; end
        end
        if md.det                          % replicate for tables / tests
            f = repmat(f, 1, cfg.runs); tm = repmat(tm, 1, cfg.runs);
            fe = repmat(fe, 1, cfg.runs); curves = repmat(curves, cfg.runs, 1);
        end
        M(k).name = md.name; M(k).role = md.role; M(k).det = md.det;
        M(k).f = f; M(k).time = tm; M(k).fe = fe;
        M(k).curve = mean(curves, 1); M(k).bestx = bestx;
        fprintf('%-10s min %7.2f  mean %8.3f  max %7.2f  std %6.3f  time %7.3fs\n', ...
            md.name, min(f), mean(f), max(f), std(f), mean(tm));
    end
    res.scen = sc; res.P = P; res.methods = M; res.grid = grid;
    res.cfg = cfg;
    save('-v7', fullfile(cfg.outdir, sprintf('raw_%s.mat', sc.name)), 'res');
end
end

function c = resample_curve(hfe, hf, grid)
% best-so-far value at every grid point (step function)
c = zeros(1, numel(grid));
for g = 1:numel(grid)
    k = find(hfe <= grid(g), 1, 'last');
    if isempty(k), k = 1; end
    c(g) = hf(k);
end
end
