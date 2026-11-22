function run_sensitivity(which)
%RUN_SENSITIVITY  One-factor-at-a-time sensitivity analysis of LD-EPO
%   parameters N, p_LA, n_LA and kappa.
%   All other parameters are kept at their default values (cfg.sens.base).
%
%   run_sensitivity                 % all parameters
%   run_sensitivity({'N','p_LA'})   % only some (to run in parallel)
%
%   Output: results/sensitivity_<param>.csv / .mat (merged by make_report)
cfg = config_experiment();
if ~exist(cfg.outdir, 'dir'), mkdir(cfg.outdir); end
sc = cfg.scen(cfg.sens.scen);
P = get_instance(cfg, cfg.sens.scen);
names = fieldnames(cfg.sens.grid);
if nargin >= 1, names = names(ismember(names, which)); end
for a = 1:numel(names)
    S = struct('param', {}, 'value', {}, 'f', {}, 'time', {});
    fid = fopen(fullfile(cfg.outdir, sprintf('sensitivity_%s.csv', names{a})), 'w');
    fprintf(fid, 'scenario,parameter,value,min,mean,max,std,time_s\n');
    for v = cfg.sens.grid.(names{a})
        prm = cfg.sens.base; prm.(names{a}) = v;
        f = zeros(1, cfg.sens.runs); tm = f;
        for r = 1:cfg.sens.runs
            rng(r, 'twister'); t0 = tic;
            o = alg_ldepo(P, prm, cfg.maxFE);
            tm(r) = toc(t0); f(r) = o.f;
        end
        S(end+1) = struct('param', names{a}, 'value', v, 'f', f, 'time', tm); %#ok<AGROW>
        fprintf('%s  %-6s = %-5g  mean %.3f  std %.3f  min %g\n', sc.name, names{a}, v, mean(f), std(f), min(f));
        fprintf(fid, '%s,%s,%g,%g,%.4f,%g,%.4f,%.3f\n', sc.name, names{a}, v, min(f), mean(f), max(f), std(f), mean(tm));
    end
    fclose(fid);
    save('-v7', fullfile(cfg.outdir, sprintf('sensitivity_%s.mat', names{a})), 'S');
end
end
