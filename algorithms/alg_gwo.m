function res = alg_gwo(P, prm, maxFE)
%ALG_GWO  Grey Wolf Optimizer (Mirjalili et al., 2014); a decreases
%   linearly 2 -> 0.  Rounding decoder.
N = getdef(prm, 'N', 30);
init = getdef(prm, 'init', 'heuristic');
[n, m] = size(P);
lb = 0.5; ub = m + 0.5 - 1e-9;
Tmax = max(1, floor((maxFE - N) / N));

Z = double(init_population(P, N, init));   % integer positions decode exactly
f = eval_pop(decode(Z, m), P); fe = N;
[fs, o] = sort(f);
Al = Z(o(1), :); Be = Z(o(min(2, N)), :); De = Z(o(min(3, N)), :);
fA = fs(1); fB = fs(min(2, N)); fD = fs(min(3, N));
hfe = fe; hf = fA;
for t = 1:Tmax
    a = 2 - 2 * t / Tmax;
    X1 = upd(Al, Z, a, N, n); X2 = upd(Be, Z, a, N, n); X3 = upd(De, Z, a, N, n);
    Z = min(max((X1 + X2 + X3) / 3, lb), ub);
    f = eval_pop(decode(Z, m), P); fe = fe + N;
    for i = 1:N                                   % update leaders
        if f(i) < fA
            fD = fB; De = Be; fB = fA; Be = Al; fA = f(i); Al = Z(i, :);
        elseif f(i) < fB
            fD = fB; De = Be; fB = f(i); Be = Z(i, :);
        elseif f(i) < fD
            fD = f(i); De = Z(i, :);
        end
    end
    hfe(end+1) = fe; hf(end+1) = fA; %#ok<AGROW>
end
res.x = decode(Al, m); res.f = fA; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function Xn = upd(L, Z, a, N, n)
A = 2 * a * rand(N, n) - a; C = 2 * rand(N, n);
D = abs(C .* repmat(L, N, 1) - Z);
Xn = repmat(L, N, 1) - A .* D;
end

function X = decode(Z, m)
X = min(max(round(Z), 1), m);
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
