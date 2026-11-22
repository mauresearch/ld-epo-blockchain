function res = alg_random_search(P, prm, maxFE)
%ALG_RANDOM_SEARCH  "Random" baseline: maxFE uniformly random assignments,
%   the best one is returned (same evaluation budget as the metaheuristics).
[n, m] = size(P);
batch = getdef(prm, 'batch', 500);
fe = 0; fbest = inf; xbest = [];
hfe = []; hf = [];
while fe < maxFE
    b = min(batch, maxFE - fe);
    X = randi(m, b, n);
    f = eval_pop(X, P); fe = fe + b;
    [fb, ib] = min(f);
    if fb < fbest, fbest = fb; xbest = X(ib, :); end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end
res.x = xbest; res.f = fbest; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
