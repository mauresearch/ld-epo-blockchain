function res = alg_epo_continuous(P, prm, maxFE)
%ALG_EPO_CONTINUOUS  Ablation variant "EPO (basic)": the ORIGINAL continuous
%   Emperor Penguin Optimizer (Dhiman & Kumar, Knowledge-Based Systems 2018)
%   applied through a simple rounding decoder: position z in [0.5, m+0.5]^n
%   -> x = round(z).  No discrete operators and no load-aware move.
%
%   Original update equations (per dimension):
%     T  = 0 if R > 1 else 1,   R = rand
%     T' = T - Max_iter / (t - Max_iter)
%     P_grid = |P_ep - P|
%     A  = M (T' + P_grid) rand - T'
%     S  = ( sqrt( f exp(-t/l) - exp(-t) ) )^2
%     C  = rand
%     D  = | S P - C P_ep |
%     P  <- P_ep - A D
%   Parameters: M = 2, f in [2,3], l in [1.5,2] (as in the original paper).
N = getdef(prm, 'N', 30);
init = getdef(prm, 'init', 'heuristic');
M = getdef(prm, 'M', 2);
fpar = getdef(prm, 'f', 2.5);
lpar = getdef(prm, 'l', 1.75);
[n, m] = size(P);
lb = 0.5; ub = m + 0.5 - 1e-9;
Tmax = max(1, floor((maxFE - N) / N));

Z = double(init_population(P, N, init));   % integer positions decode exactly
X = decode(Z, m);
f = eval_pop(X, P); fe = N;
[fbest, ib] = min(f); zbest = Z(ib, :); xbest = X(ib, :);
hfe = fe; hf = fbest;
for t = 1:Tmax
    R = rand; T = double(~(R > 1));
    Tp = T - Tmax / ((t - Tmax) - 1e-9);     % avoid division by zero at t = Tmax
    Tp = min(Tp, 10);                         % numeric safeguard
    S = (sqrt(max(fpar * exp(-t / lpar) - exp(-t), 0)))^2;
    for i = 1:N
        Pg = abs(zbest - Z(i, :));
        A = M * (Tp + Pg) .* rand(1, n) - Tp;
        C = rand(1, n);
        D = abs(S * Z(i, :) - C .* zbest);
        Z(i, :) = zbest - A .* D;
    end
    Z = min(max(Z, lb), ub);
    X = decode(Z, m);
    f = eval_pop(X, P); fe = fe + N;
    [fb, ib] = min(f);
    if fb < fbest, fbest = fb; zbest = Z(ib, :); xbest = X(ib, :); end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end
res.x = xbest; res.f = fbest; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function X = decode(Z, m)
X = min(max(round(Z), 1), m);
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
