% handle-1 : run 
% varying nof kernels w/ fixed nof examples

%% baseline 1: svm w/ sum kernel
%play_svm;
play_svm_largescale;
%
%% baseline 2: 1-norm cplex
%play_lp_norm(1,1);
play_lp_norm(1,1,0);
%
%% baseline 3: simple mkl
%play_simple_mkl;
%
%% cplex, lp-norm
%play_lp_norm(1.333,1);
%play_lp_norm(2,1);
%play_lp_norm(3,1);
%
%% direct, lp-norm
play_lp_norm(1,3);
%play_lp_norm(1.333,3);
%play_lp_norm(2,3);
%play_lp_norm(3,3);
%
%% newton, lp-norm
%play_lp_norm(1.333,2);
%play_lp_norm(2,2);
%play_lp_norm(3,2);
%
%% cplex, large scale
%play_lp_norm_largescale(1,1);
%play_lp_norm_largescale(1.333,1);
%play_lp_norm_largescale(2,1);
%play_lp_norm_largescale(3,1);
%
%% direct, large scale
%play_lp_norm_largescale(1.333,3);
%play_lp_norm_largescale(2,3);
%play_lp_norm_largescale(3,3);
%
%% newton, large scale
%play_lp_norm_largescale(1.333,2);
%play_lp_norm_largescale(2,2);
%play_lp_norm_largescale(3,2);

%play_lp_norm(1,1);
%play_lp_norm_largescale(2,1);
play_lp_norm_largescale(2,3);
