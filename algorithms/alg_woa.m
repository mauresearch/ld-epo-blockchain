function res = alg_woa(P, prm, maxFE)
%ALG_WOA  Whale Optimization Algorithm (Mirjalili & Lewis, 2016),
%   spiral constant b = 1, a decreases 2 -> 0.  Rounding decoder.
N = getdef(prm, 'N', 30);
init = getdef(prm, 'init', 'heuristic');
b = getdef(prm, 'b', 1);
[n, m] = size(P);
lb = 0.5; ub = m + 0.5 - 1e-9;
Tmax = max(1, floor((maxFE - N) / N));

Z = double(init_population(P, N, init));   % integer positions decode exactly
f = eval_pop(decode(Z, m), P); fe = N;
[fbest, ib] = min(f); zb = Z(ib, :);
hfe = fe; hf = fbest;
for t = 1:Tmax
    a = 2 - 2 * t / Tmax;
    for i = 1:N
        r1 = rand; r2 = rand; A = 2 * a * r1 - a; C = 2 * r2;
        p = rand; l = 2 * rand - 1;
        if p < 0.5
            if abs(A) < 1
                D = abs(C * zb - Z(i, :));       Z(i, :) = zb - A * D;
            else
                zr = Z(randi(N), :);
                D = abs(C * zr - Z(i, :));       Z(i, :) = zr - A * D;
            end
        else
            D = abs(zb - Z(i, :));
            Z(i, :) = D * exp(b * l) * cos(2 * pi * l) + zb;
        end
    end
    Z = min(max(Z, lb), ub);
    f = eval_pop(decode(Z, m), P); fe = fe + N;
    [fb, ib] = min(f);
    if fb < fbest, fbest = fb; zb = Z(ib, :); end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end
res.x = decode(zb, m); res.f = fbest; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function X = decode(Z, m)
X = min(max(round(Z), 1), m);
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
