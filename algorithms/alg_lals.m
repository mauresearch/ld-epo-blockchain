function res = alg_lals(P, prm, maxFE)
%ALG_LALS  Ablation variant "LA-LS": population of independent local
%   searches using the SAME neighbourhood ingredients as LD-EPO (relocation
%   of one block with relocate_genes.m + load-aware move with probability
%   p_LA, greedy acceptance) but WITHOUT the EPO mechanics (no emperor /
%   huddle move, no temperature / |A| based exploration-exploitation
%   switching, no information sharing between penguins).
%   Used to show what EPO contributes beyond a plain local search.
N    = getdef(prm, 'N', 30);
p_LA = getdef(prm, 'p_LA', 1.0);
n_LA = getdef(prm, 'n_LA', 3);
init = getdef(prm, 'init', 'heuristic');
[n, m] = size(P);

X = init_population(P, N, init);
[f, L] = eval_pop(X, P);
fe = N;
[fbest, ib] = min(f); xbest = X(ib, :);
hfe = fe; hf = fbest;
while fe < maxFE
    for i = 1:N
        [y, Ly] = relocate_genes(X(i, :), L(i, :), P, 1);  fe = fe + 1;
        if rand < p_LA
            for r = 1:n_LA
                [y, Ly, mv] = load_aware_move(y, Ly, P);
                if ~mv, break; end
                fe = fe + 1;
            end
        end
        fy = max(Ly);
        if fy <= f(i)
            X(i, :) = y; f(i) = fy; L(i, :) = Ly;
            if fy < fbest, fbest = fy; xbest = y; end
        end
    end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end
res.x = xbest; res.f = fbest; res.fe = fe; res.hist_fe = hfe; res.hist_f = hf;
res.info = struct();
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
