function res = alg_ga(P, prm, maxFE)
%ALG_GA  Integer-coded genetic algorithm.
%   Encoding: same assignment vector as LD-EPO.
%   Binary tournament selection, uniform crossover (pc), per-gene random
%   reset mutation (pm = 1/n by default), elitism (n_elite).
N  = getdef(prm, 'N', 30);
pc = getdef(prm, 'pc', 0.9);
pm = getdef(prm, 'pm', []);
ne = getdef(prm, 'n_elite', 2);
init = getdef(prm, 'init', 'heuristic');
[n, m] = size(P);
if isempty(pm), pm = 1 / n; end

X = init_population(P, N, init);
f = eval_pop(X, P); fe = N;
[fbest, ib] = min(f); xbest = X(ib, :);
hfe = fe; hf = fbest;
while fe < maxFE
    nc = min(N - ne, maxFE - fe);
    [~, ord] = sort(f);
    E = X(ord(1:ne), :); fE = f(ord(1:ne));
    Cc = zeros(nc, n);
    for k = 1:2:nc
        p1 = tour(f); p2 = tour(f);
        a = X(p1, :); b = X(p2, :);
        if rand < pc
            mask = rand(1, n) < 0.5;
            c1 = a; c1(mask) = b(mask);
            c2 = b; c2(mask) = a(mask);
        else
            c1 = a; c2 = b;
        end
        Cc(k, :) = c1;
        if k + 1 <= nc, Cc(k + 1, :) = c2; end
    end
    mut = rand(nc, n) < pm;
    R = randi(m, nc, n);
    Cc(mut) = R(mut);
    fc = eval_pop(Cc, P); fe = fe + nc;
    X = [E; Cc]; f = [fE; fc];
    if size(X, 1) < N                     % last partial generation
        X = [X; repmat(xbest, N - size(X, 1), 1)];
        f = [f; repmat(fbest, N - numel(f), 1)];
    end
    [fb, ib] = min(f);
    if fb < fbest, fbest = fb; xbest = X(ib, :); end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end
res.x = xbest; res.f = fbest; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function w = tour(ff)
c = randi(numel(ff), 1, 2);
if ff(c(1)) <= ff(c(2)), w = c(1); else, w = c(2); end
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
