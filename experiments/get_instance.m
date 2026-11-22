function P = get_instance(cfg, s)
%GET_INSTANCE  Processing-time matrix of scenario s.
%   If instances/P_<name>.csv exists (shipped with the code), it is read,
%   so that MATLAB and Octave use exactly the same test instances (their
%   random generators differ).  Otherwise generate_instance is called.
sc = cfg.scen(s);
here = fileparts(fileparts(mfilename('fullpath')));   % repository root
fn = fullfile(here, 'instances', sprintf('P_%s.csv', sc.name));
if exist(fn, 'file')
    P = dlmread(fn);
    if ~isequal(size(P), [sc.n, sc.m])
        error('instance file %s has wrong size', fn);
    end
else
    inst = generate_instance(sc.n, sc.m, sc.seed, cfg.gen);
    P = inst.P;
end
end
