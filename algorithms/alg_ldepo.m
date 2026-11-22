function res = alg_ldepo(P, prm, maxFE)
%ALG_LDEPO  Discrete, Load-aware Emperor Penguin Optimizer (LD-EPO).
%
%   res = alg_ldepo(P, prm, maxFE)
%
%   prm fields (defaults in brackets):
%     N     population size (number of penguins)            [30]
%     p_LA  probability of applying the load-aware move     [1.0]
%     n_LA  max. number of load-aware moves per application [3]
%     M     EPO movement parameter                          [2]
%     kappa exploration step: genes relocated = ceil(kappa |A|) [2]
%     init  'heuristic' (LPT-MCT + Min-Min + random) or 'random' ['heuristic']
%   With p_LA = 0 the algorithm is the ablation variant "D-EPO"
%   (discrete EPO without the load-aware operator).
%
%   Mapping of the continuous EPO (Dhiman & Kumar, 2018) to the discrete
%   assignment space:
%
%   continuous EPO                      discrete LD-EPO
%   -------------------------------     --------------------------------------
%   position P(x) in R^n                assignment vector X in {1..m}^n
%   emperor (best) P_ep                 X_best
%   P_grid = |P - P_ep|                 G_i = Hamming(X_i, X_best) / n
%   temperature T' (huddle)             T'_t = 1 - FE_t / FE_max (cooling huddle)
%   A = M (T' + P_grid) rand - T'       A = M (T'_t + G_i) r1 - T'_t
%   C = rand                            C = r2
%   D_ep = |S(A) P - C P_ep|            D_i = { k : X_i(k) ~= X_best(k) }
%   P <- P_ep - A D_ep                  if |A| <  1 : huddle (exploitation):
%                                           q = ceil(C (1-|A|) |D_i|) blocks of
%                                           D_i move to the emperor's miner if
%                                           the makespan is not exceeded
%                                           (huddle_move.m)
%                                       if |A| >= 1 : relocation (exploration):
%                                           relocate k = ceil(kappa |A|) random
%                                           blocks to miners that can take them
%                                           without exceeding the makespan
%                                           (relocate_genes.m)
%   (no explicit operator)              load-aware move with prob. p_LA
%   (no selection)                      non-worsening (greedy) acceptance
%
%   Early in the search T'_t ~ 1 and G_i is large, so |A| >= 1 is frequent
%   (exploration).  As t -> T_max and the huddle contracts (G_i -> 0),
%   |A| -> 0 and the penguins exploit the neighbourhood of the emperor.
%   This is what distinguishes LD-EPO from a plain random local search
%   (variant "LA-LS" in the ablation study).

N    = getdef(prm, 'N', 30);
p_LA = getdef(prm, 'p_LA', 1.0);
n_LA = getdef(prm, 'n_LA', 3);
M    = getdef(prm, 'M', 2);
kappa = getdef(prm, 'kappa', 2);
init = getdef(prm, 'init', 'heuristic');

[n, m] = size(P);

X = init_population(P, N, init);
[f, L] = eval_pop(X, P);
fe = N;
[fbest, ib] = min(f); xbest = X(ib, :);
hfe = fe; hf = fbest;
n_explore = 0; n_exploit = 0; n_la_moved = 0; n_la_tried = 0;

while fe < maxFE
    Tp = max(0, 1 - fe / maxFE);               % huddle temperature T' (budget-based t/T_max)
    for i = 1:N
        xi = X(i, :);
        D = find(xi ~= xbest);
        G = numel(D) / n;                      % P_grid analogue
        A = M * (Tp + G) * rand - Tp;
        C = rand;
        if abs(A) >= 1 || isempty(D)           % ---- exploration (relocation)
            k = ceil(kappa * max(abs(A), 1));
            [y, Ly] = relocate_genes(xi, L(i, :), P, k);
            n_explore = n_explore + 1;
        else                                   % ---- exploitation (huddle)
            q = ceil(C * (1 - abs(A)) * numel(D));
            pos = D(randperm(numel(D), q));
            [y, Ly] = huddle_move(xi, L(i, :), xbest, pos, P);
            n_exploit = n_exploit + 1;
        end
        fe = fe + 1;
        if rand < p_LA                          % ---- load-aware move
            for r = 1:n_LA
                n_la_tried = n_la_tried + 1;
                [y, Ly, mv] = load_aware_move(y, Ly, P);
                if ~mv, break; end
                n_la_moved = n_la_moved + 1;
                fe = fe + 1;                    % counted as one evaluation
            end
        end
        fy = max(Ly);
        if fy <= f(i)                           % greedy (non-worsening) acceptance
            X(i, :) = y; f(i) = fy; L(i, :) = Ly;
            if fy < fbest
                fbest = fy; xbest = y;
            end
        end
    end
    hfe(end+1) = fe; hf(end+1) = fbest; %#ok<AGROW>
end

res.x = xbest; res.f = fbest; res.fe = fe;
res.hist_fe = hfe; res.hist_f = hf;
res.info = struct('n_explore', n_explore, 'n_exploit', n_exploit, ...
                  'la_success_rate', n_la_moved / max(1, n_la_tried));
end

function v = getdef(s, f, d)
if isfield(s, f), v = s.(f); else, v = d; end
end
