num=16;
i=0;

clf



clear mkl;
i=i+1;
load ../nofex/results/res_cplex_1.000norm_wrapper.mat
prep_plot(mkl.timex, 'm-', i, num, 'm'); % 1-norm wrapper

clear smkl;
i=i+1;
load ../nofex/results/res_simplemkl.mat
prep_plot(smkl.timex, 'm--', i, num, 'm'); % 1-norm simpleMKL 



clear mkl_ls;
i=i+1;
load ../nofex/results/res_cplex_1.000norm_ls.mat
prep_plot(mkl_ls.timex, 'r-', i, num) %MKL 1-norm SIP

clear mkl_ls;
i=i+1;
load ../nofex/results/res_cplex_2.000norm_ls.mat
prep_plot(mkl_ls.timex, 'r:', i, num) %MKL 2-norm SIP

clear mkl_ls;
i=i+1;
load ../nofex/results/res_direct_2.000norm_ls.mat
prep_plot(mkl_ls.timex, 'r--', i, num) %MKL 2-norm direct

clear svm_ls;
i=i+1;
load ../nofex/results/res_svm_ls.mat
prep_plot(svm_ls.timex, 'r.', i, num) %SVM



clear hmkl;
i=i+1;
load ../nofex/results/res_hessmkl.mat
prep_plot(hmkl.timex, 'r-.', i, num); % 1-norm hessMKL



clear mkl;
i=i+1;
load ../nofex/results/res_direct_1.000norm.mat
prep_plot(mkl.timex, 'b-', i, num);% 1.000-norm direct

clear mkl;
i=i+1;
load ../nofex/results/res_direct_1.333norm.mat
prep_plot(mkl.timex, 'b--', i, num);% 1.333-norm direct

clear mkl;
i=i+1;
load ../nofex/results/res_direct_2.000norm.mat
prep_plot(mkl.timex, 'b:', i, num); % 2-norm direct

clear mkl;
i=i+1;
load ../nofex/results/res_direct_3.000norm.mat
prep_plot(mkl.timex, 'b-.', i, num); % 3-norm direct



clear mkl;
i=i+1;
load ../nofex/results/res_cplex_1.000norm.mat % 1-norm cplex
prep_plot(mkl.timex, 'r-', i, num);

clear mkl;
i=i+1;
load ../nofex/results/res_cplex_1.333norm.mat
prep_plot(mkl.timex, 'r--', i, num);% 1.333-norm cplex

clear mkl;
i=i+1;
load ../nofex/results/res_cplex_2.000norm.mat
prep_plot(mkl.timex, 'r:', i, num); % 2-norm cplex

clear mkl;
i=i+1;
load ../nofex/results/res_cplex_3.000norm.mat
prep_plot(mkl.timex, 'r-.', i, num); % 3-norm cplex

%clear mkl;
%load ../nofex/results/res_newton_1.333norm.mat
%prep_plot(mkl.timex, 'g--', i, num);% 1.333-norm newton
%
%clear mkl;
%load ../nofex/results/res_newton_2.000norm.mat
%prep_plot(mkl.timex, 'g:', i, num); % 2-norm newton
%
%clear mkl;
%load ../nofex/results/res_newton_3.000norm.mat
%prep_plot(mkl.timex, 'g-.', i, num); % 3-norm newton

clear svm;
i=i+1;
load ../nofex/results/res_svm.mat
prep_plot(svm.timex, 'k-', i, num) %SVM



title('50 precomputed kernels - loglog plot')
ylabel('training time (seconds)')
xlabel('number of training examples')
legend( ...
	 'MKL wrapper  (1-norm SILP)', 'SimpleMKL', ...
		'on-the-fly SILP       (1-norm)', '               SIP         (2-norm)', '               analytical (2-norm)',  ...
		'               SVM        (\infty-norm)', ...
		'HessianMKL',...
		'interleaved analytical (1-norm)', '                                 (4/3-norm)','                                 (2-norm)','                                 (3-norm)', ...
		'interleaved SIP (1-norm)', '                         (4/3-norm)','                         (2-norm)','                         (3-norm)', ...
		            'SVM  (\infty-norm)',...
	 'Location','SouthEastOutside');
%'1.333-norm newton','2-norm newton','3-norm newton', ...
%legend('SVM','1-norm SILP','1-norm simpleMKL','2-norm Newton','3-norm Newton','2-norm cutting plane','3-norm cutting plane');

set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
grid on;
axis tight;
