#!/usr/bin/env python3
"""Exact MILP for the block-to-miner assignment problem (R||Cmax), solved with HiGHS.

    min C  s.t.  sum_j x_ij = 1  (all i),  sum_i p_ij x_ij <= C  (all j),  x_ij in {0,1}

Usage:
    pip install highspy numpy
    python tools/milp_highs.py instances/P_S1.csv [time_limit_s]
    python tools/milp_highs.py --all [time_limit_s]      # writes results/milp_results.csv

This script produced results/milp_results.csv used in the paper (HiGHS 1.16,
600 s, single thread). In MATLAB, run_milp.m solves the same model with intlinprog.
"""
import sys, time, glob, os
import numpy as np
import highspy


def solve(P, time_limit=600.0):
    n, m = P.shape
    h = highspy.Highs()
    h.setOptionValue('output_flag', False)
    h.setOptionValue('threads', 1)
    h.setOptionValue('time_limit', float(time_limit))
    inf = highspy.kHighsInf
    nv = n * m + 1                      # x(i,j) -> j*n + i ; C -> last
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = nv, n + m
    lp.col_cost_ = np.r_[np.zeros(n * m), 1.0]
    lp.col_lower_ = np.zeros(nv)
    lp.col_upper_ = np.r_[np.ones(n * m), inf]
    lp.row_lower_ = np.r_[np.ones(n), -inf * np.ones(m)]
    lp.row_upper_ = np.r_[np.ones(n), np.zeros(m)]
    starts, idx, val = [], [], []
    for j in range(m):
        for i in range(n):
            starts.append(len(idx)); idx += [i, n + j]; val += [1.0, float(P[i, j])]
    starts.append(len(idx)); idx += list(range(n, n + m)); val += [-1.0] * m
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = starts + [len(idx)]
    lp.a_matrix_.index_, lp.a_matrix_.value_ = idx, val
    lp.integrality_ = [highspy.HighsVarType.kInteger] * (n * m) + [highspy.HighsVarType.kContinuous]
    h.passModel(lp)
    t0 = time.time(); h.run(); t = time.time() - t0
    info = h.getInfo()
    status = 'optimal' if h.getModelStatus() == highspy.HighsModelStatus.kOptimal else 'time_limit'
    return info.objective_function_value, info.mip_dual_bound, status, t


if __name__ == '__main__':
    args = sys.argv[1:]
    if not args:
        print(__doc__); sys.exit(1)
    if args[0] == '--all':
        tl = float(args[1]) if len(args) > 1 else 600
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        rows = ['scenario,m,n,milp_best,milp_lower_bound,status,time_s']
        for f in sorted(glob.glob(os.path.join(root, 'instances', 'P_S*.csv'))):
            P = np.loadtxt(f, delimiter=','); n, m = P.shape
            name = os.path.basename(f)[2:-4]
            best, lb, st, t = solve(P, tl)
            print(f'{name}: best={best:g} lb={lb:g} {st} {t:.1f}s')
            rows.append(f'{name},{m},{n},{best:g},{lb:g},{st},{t:.2f}')
        os.makedirs(os.path.join(root, 'results'), exist_ok=True)
        open(os.path.join(root, 'results', 'milp_results.csv'), 'w').write('\n'.join(rows) + '\n')
    else:
        P = np.loadtxt(args[0], delimiter=',')
        tl = float(args[1]) if len(args) > 1 else 600
        best, lb, st, t = solve(P, tl)
        print(f'{args[0]}: best={best:g} lower_bound={lb:g} status={st} time={t:.1f}s')
