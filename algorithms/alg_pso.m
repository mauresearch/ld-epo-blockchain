function res = alg_pso(P, prm, maxFE)
%ALG_PSO  Particle swarm optimisation (Kennedy & Eberhart 1995) with
%   linearly decreasing inertia (0.9 -> 0.4), c1 = c2 = 2, velocity clamp
%   vmax = 0.2 (ub - lb).  Continuous positions in [0.5, m+0.5]^n are
%   decoded by rounding (same decoder as all continuous baselines).
N = getdef(prm, 'N', 30);
init = getdef(prm, 'init', 'heuristic');
c1 = getdef(prm, 'c1', 2); c2 = getdef(prm, 'c2', 2);
w1 = getdef(prm, 'w_max', 0.9); w2 = getdef(prm, 'w_min', 0.4);
[n, m] = size(P);
lb = 0.5; ub = m + 0.5 - 1e-9; vmax = 0.2 * (ub - lb);
Tmax = max(1, floor((maxFE - N) / N));

Z = double(init_population(P, N, init)); V = zeros(N, n);   % integer positions decode exactly
f = eval_pop(decode(Z, m), P); fe = N;
Pb = Z; fPb = f;
[fbest, ib] = min(f); gz = Z(ib, :);
hfe = fe; hf = fbest;
for t = 1:Tmax
    w = w1 - (w1 - w2) * t / Tmax;
    V = w * V + c1 * rand(N, n) .* (Pb - Z) + c2 * rand(N, n) .* bsxfun(@minus, gz, Z);
    V = min(max(V, -vmax), vmax);
    Z = min(max(Z + V, lb), ub);
    f = eval_pop(decode(Z, m), P); fe = fe + N;
    imp = f < fPb; Pb(imp, :) = Z(imp, :); fPb(imp) = f(imp);
    [fb, ib] = min(fPb);
    if fb < fbest, fbest = fb; gz = Pb(ib, :); end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end
res.x = decode(gz, m); res.f = fbest; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function X = decode(Z, m)
X = min(max(round(Z), 1), m);
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
