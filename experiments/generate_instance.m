function inst = generate_instance(n, m, seed, opts)
%GENERATE_INSTANCE  Reproducible block-to-miner instance (processing-time matrix P).
%
%   inst = generate_instance(n, m, seed)
%   inst = generate_instance(n, m, seed, opts)
%
%   Model:
%     s_i   : workload of block i (e.g. number/size of transactions to verify),
%             integer ~ U{s_min, ..., s_max}
%     v_j   : relative computing capacity of miner j ~ U(v_min, v_max)
%             (miner heterogeneity)
%     e_ij  : block-miner affinity noise ~ U(1-eps, 1+eps) that models
%             bandwidth / verification / hardware effects that are specific
%             to the (block, miner) pair -> unrelated-machine model R||Cmax
%     p_ij  = max(1, round( s_i / v_j * e_ij ))
%
%   Default parameters: s in [5,20], v in [0.5,2.0], eps = 0.10.
%   The same (n, m, seed) always yields the same matrix.

if nargin < 4, opts = struct(); end
s_min = getdef(opts, 's_min', 5);
s_max = getdef(opts, 's_max', 20);
v_min = getdef(opts, 'v_min', 0.5);
v_max = getdef(opts, 'v_max', 2.0);
eps_  = getdef(opts, 'eps', 0.10);

old = rng;               % do not disturb caller's random stream
rng(seed, 'twister');
s = randi([s_min, s_max], n, 1);
v = v_min + (v_max - v_min) * rand(1, m);
e = 1 - eps_ + 2 * eps_ * rand(n, m);
P = max(1, round(bsxfun(@rdivide, s, v) .* e));
rng(old);

inst.n = n;  inst.m = m;  inst.seed = seed;
inst.P = P;  inst.s = s;  inst.v = v;
inst.params = struct('s_min', s_min, 's_max', s_max, 'v_min', v_min, ...
                     'v_max', v_max, 'eps', eps_);
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
