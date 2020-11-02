num=16;
i=0;

clf

clear mkl;
i=i+1;
load ../nofks/results/res_cplex_1.000norm_wrapper.mat
prep_plot(mkl.timex, 'm:', i, num); % 1-norm wrapper

clear smkl;
i=i+1;
load ../nofks/results/res_simplemkl.mat
prep_plot(smkl.timex, 'm-', i, num); % 1-norm simpleMKL 




i=i+1;

clear mkl_ls;
i=i+1;
load ../nofks/results/res_cplex_2.000norm_ls.mat
prep_plot(mkl_ls.timex, 'c-', i, num); % 2-norm cplex ls

clear mkl_ls;
i=i+1;
load ../nofks/results/res_direct_2.000norm_ls.mat
prep_plot(mkl_ls.timex, 'c:', i, num); % 2-norm direct ls

clear svm_ls;
i=i+1;
load ../nofks/results/res_svm_ls.mat
prep_plot(svm_ls.timex, 'k:', i, num) %SVM





%clear mkl;
%load ../nofks/results/res_newton_1.333norm.mat
%prep_plot(mkl.timex, 'g--');% 1.333-norm newton
%
%clear mkl;
%load ../nofks/results/res_newton_2.000norm.mat
%prep_plot(mkl.timex, 'g:'); % 2-norm newton
%
%clear mkl;
%load ../nofks/results/res_newton_3.000norm.mat
%prep_plot(mkl.timex, 'g-.'); % 3-norm newton
%
%



clear hmkl;
i=i+1;
load ../nofks/results/res_hessmkl.mat
prep_plot(hmkl.timex, 'm-', i, num); % 1-norm hessMKL 





clear mkl;
i=i+1;
load ../nofks/results/res_direct_1.000norm.mat
prep_plot(mkl.timex, 'b-', i, num);% 1.000-norm direct

clear mkl;
i=i+1;
load ../nofks/results/res_direct_1.333norm.mat
prep_plot(mkl.timex, 'b--', i, num);% 1.333-norm direct

clear mkl;
i=i+1;
load ../nofks/results/res_direct_2.000norm.mat
prep_plot(mkl.timex, 'b:', i, num); % 2-norm direct

clear mkl;
i=i+1;
load ../nofks/results/res_direct_3.000norm.mat
prep_plot(mkl.timex, 'b-.', i, num); % 3-norm direct






clear mkl;
i=i+1;
load ../nofks/results/res_cplex_1.000norm.mat % 1-norm cplex
prep_plot(mkl.timex, 'r-', i, num);

clear mkl;
i=i+1;
load ../nofks/results/res_cplex_1.333norm.mat
prep_plot(mkl.timex, 'r--', i, num);% 1.333-norm cplex

clear mkl;
i=i+1;
load ../nofks/results/res_cplex_2.000norm.mat
prep_plot(mkl.timex, 'r:', i, num); % 2-norm cplex

clear mkl;
i=i+1;
load ../nofks/results/res_cplex_3.000norm.mat
prep_plot(mkl.timex, 'r-.', i, num); % 3-norm cplex



clear svm;
i=i+1;
load ../nofks/results/res_svm.mat
prep_plot(svm.timex, 'k-', i, num) %SVM



title('1000 examples with varying number of precomputed kernels')
ylabel('training time (seconds; logarithmic)')
xlabel('number of kernels')
%legend('SVM','SVM-LS', '1-norm simpleMKL', ...
%		'1-norm SILP', '1.333-norm SIP','2-norm SIP','3-norm SIP', ...
%		'1.333-norm newton','2-norm newton','3-norm newton', ...
%		'1.333-norm direct','2-norm direct','3-norm direct', 'Location','NorthEastOutside');
legend( ...
	 'MKL wrapper  (1-norm SILP)', 'SimpleMKL', ...
		'on-the-fly SIP         (2-norm)', '               analytical (2-norm)',  ...
		'               SVM        (\infty-norm)', ...
		'HessianMKL',...
		'interleaved analytical (1-norm)', '                                 (4/3-norm)','                                 (2-norm)','                                 (3-norm)', ...
		'interleaved SIP (1-norm)', '                         (4/3-norm)','                         (2-norm)','                         (3-norm)', ...
		            'SVM  (\infty-norm)',...
	 'Location','NorthEastOutside');
%legend('1-norm SILP wrapper', '1-norm simpleMKL', '1-norm hessianMKL',...
%		'2-norm SIP on-the-fly', '2-norm direct on-the-fly', 'SVM on-the-fly', ...
%		'3-norm SILP', '2-norm SIP', '1.333-norm SIP', '1-norm SILP', ...
%		'3-norm direct', '2-norm direct','1.333-norm direct','1-norm direct',...
%		'SVM', ...
%		'Location','SouthEastOutside');%,'2-norm Newton','3-norm Newton','2-norm cutting plane','3-norm cutting plane');

set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
grid on;
axis tight;
