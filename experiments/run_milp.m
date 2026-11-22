function run_milp()
%RUN_MILP  Exact MILP solution of the scenarios (optimality gap).
%   Writes results/milp_results.csv : scenario, best integer value,
%   proven lower bound, status and CPU time.  The best LD-EPO solution of
%   run_experiments is passed as a starting incumbent when available.
cfg = config_experiment();
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end
fid = fopen(fullfile(cfg.outdir, 'milp_results.csv'), 'w');
fprintf(fid, 'scenario,m,n,milp_best,milp_lower_bound,status,time_s\n');
for s = cfg.milp_scen
    sc = cfg.scen(s);
    P = get_instance(cfg, s);
    x0 = [];
    fr = fullfile(cfg.outdir, sprintf('raw_%s.mat', sc.name));
    if exist(fr, 'file')
        S = load(fr); k = find(strcmp({S.res.methods.name}, 'LD-EPO'));
        x0 = S.res.methods(k).bestx;
    end
    o = solve_milp(P, cfg.milp_time, x0);
    fprintf('%s: MILP best = %g, LB = %g (%s, %.1fs)\n', sc.name, o.f, o.lb, o.status, o.time);
    fprintf(fid, '%s,%d,%d,%g,%g,%s,%.2f\n', sc.name, sc.m, sc.n, o.f, o.lb, o.status, o.time);
end
fclose(fid);
end
