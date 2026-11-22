function X = init_population(P, N, mode)
%INIT_POPULATION  Initial population of N assignment vectors.
%   mode = 'random'    : every gene uniform in {1..m}
%   mode = 'heuristic' : row 1 = LPT-MCT solution, row 2 = Min-Min solution,
%                        the remaining N-2 rows are uniform random.
%   The SAME function is used by every population-based method so that all
%   metaheuristics start from identical information (fair comparison).
[n, m] = size(P);
X = randi(m, N, n);
if strcmpi(mode, 'heuristic')
    r1 = alg_greedy_mct(P, struct(), 1);
    X(1, :) = r1.x;
    if N >= 2
        r2 = alg_greedy_minmin(P, struct(), 1);
        X(2, :) = r2.x;
    end
end
end
