function out = solve_milp(P, timelimit, x0)
%SOLVE_MILP  Exact MILP model of the block-to-miner problem 
%
%   min  C
%   s.t. sum_j x_ij = 1                  for all i      (Eq. 5)
%        sum_i p_ij x_ij - C <= 0        for all j      (Eqs. 2-3)
%        x_ij in {0,1}, C >= 0                          (Eq. 6)
%
%   Uses intlinprog (MATLAB Optimization Toolbox) if available, otherwise
%   GLPK (built into GNU Octave).
%   out.f      : best integer objective found
%   out.lb     : proven lower bound (= out.f when solved to optimality)
%   out.status : 'optimal' | 'time_limit' | 'failed'
%   out.time   : CPU seconds
%   x0 (optional) : a feasible assignment used as an upper bound / start.
if nargin < 2 || isempty(timelimit), timelimit = 600; end
if nargin < 3, x0 = []; end
[n, m] = size(P);
nv = n * m + 1;                          % x(i,j) -> (j-1)*n + i,  C -> nv
c = [zeros(n * m, 1); 1];
Aeq = [repmat(speye(n), 1, m), sparse(n, 1)];
beq = ones(n, 1);
Ain = sparse(m, nv);
for j = 1:m
    Ain(j, (j - 1) * n + (1:n)) = P(:, j)';  %#ok<SPRIX>
end
Ain(:, nv) = -1;
bin = zeros(m, 1);
lb = zeros(nv, 1); ub = [ones(n * m, 1); inf];

% simple valid lower bound used when the solver gives none
LB0 = max([max(min(P, [], 2)), ceil(sum(min(P, [], 2)) / m)]);

t0 = tic;
out = struct('f', inf, 'lb', LB0, 'status', 'failed', 'time', 0, 'x', []);
if exist('intlinprog', 'file') == 2
    opts = optimoptions('intlinprog', 'Display', 'off', 'MaxTime', timelimit);
    if ~isempty(x0)
        z0 = zeros(nv, 1);
        z0((x0(:) - 1) * n + (1:n)') = 1;
        z0(nv) = max(accumarray(x0(:), P(sub2ind(size(P), (1:n)', x0(:))), [m 1]));
        args = {z0};
    else
        args = {};
    end
    [z, fval, flag, outp] = intlinprog(c, 1:n*m, Ain, bin, Aeq, beq, lb, ub, args{:}, opts);
    if ~isempty(z)
        out.f = fval;
        out.x = decode(z, n, m);
    end
    if isfield(outp, 'relativegap') && ~isempty(z)
        out.lb = max(LB0, fval * (1 - outp.relativegap));
    end
    if flag == 1, out.status = 'optimal'; out.lb = fval;
    elseif ~isempty(z), out.status = 'time_limit'; end
else
    A = [Aeq; Ain]; b = [beq; bin];
    ctype = [repmat('S', n, 1); repmat('U', m, 1)]';
    vtype = [repmat('I', n * m, 1); 'C']';
    param.msglev = 0;
    param.tmlim = round(timelimit * 1000);
    param.presol = 1;
    [z, fmin, errnum, extra] = glpk(c, A, b, lb, ub, ctype, vtype, 1, param); %#ok<ASGLU>
    if ~isempty(z) && all(isfinite(z)) && fmin < inf && extra.status ~= 1
        out.f = fmin; out.x = decode(z, n, m);
    end
    if extra.status == 5
        out.status = 'optimal'; out.lb = fmin;
    elseif isfinite(out.f)
        out.status = 'time_limit';
    end
end
out.time = toc(t0);
if ~isempty(x0)
    fx0 = eval_pop(x0(:)', P);
    if fx0 < out.f, out.f = fx0; out.x = x0(:)'; end
    if strcmp(out.status, 'failed'), out.status = 'time_limit'; end
end
end

function x = decode(z, n, m)
X = reshape(round(z(1:n*m)), n, m);
[~, x] = max(X, [], 2); x = x';
end
