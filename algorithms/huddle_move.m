function [y, Ly] = huddle_move(y, Ly, xbest, pos, P)
%HUDDLE_MOVE  Exploitation (huddling) operator of LD-EPO.
%   The penguin moves towards the emperor X_best: for every block b in
%   'pos' (a random subset of the blocks on which y and X_best differ),
%   b is re-assigned to the emperor's miner X_best(b) provided this does
%   not push that miner above the current makespan of y.  The penguin
%   therefore inherits the emperor's structure without losing its own
%   solution quality (makespan-preserving path towards the emperor).
Cmax = max(Ly);
for b = pos
    cur = y(b); j = xbest(b);
    if Ly(j) + P(b, j) <= Cmax
        Ly(cur) = Ly(cur) - P(b, cur);
        Ly(j) = Ly(j) + P(b, j);
        y(b) = j;
    end
end
end
