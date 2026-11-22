% MAIN  Reproduces every result of the paper.
%
%   1) run_experiments : 12 methods x 6 scenarios x 30 independent runs
%   2) run_milp        : exact MILP model of every scenario (optimality gap)
%   3) run_sensitivity : parameter sensitivity of LD-EPO on S3
%   4) make_report     : all tables (results/*.csv) and figures (figures/*.png)
%
%   Quick smoke test (a few minutes):  setenv('EPO_QUICK','1'); main
%   Settings: config_experiment.m
%
%   The scenarios are independent: on a multi-core PC you can open several
%   MATLAB sessions and run e.g. run_experiments([1 2 3]) in one and
%   run_experiments([4 5 6]) in another, then run_milp / run_sensitivity /
%   make_report.
clc;
root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
cd(root);
t0 = tic;
run_experiments();
run_milp();
run_sensitivity();
make_report();
fprintf('Total time: %.1f min\n', toc(t0) / 60);
