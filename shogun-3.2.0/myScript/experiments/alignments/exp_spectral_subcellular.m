tic;
load('data/subcellular.mat');
time_data=toc;

K = kscale(K);
%y = sign(y-1.5);
yaux = y(:); 
y = zeros(size(yaux'));
y(1,find(yaux==0)) = 1;
y(2,find(yaux==1)) = 1;
y(3,find(yaux==2)) = 1;
y(4,find(yaux==3)) = 1;
y = y';

alignments = compute_alignments(K);
specs = compute_specs(alignments);
%talignments = compute_mc_targetalignments(K,y);
talignments = compute_targetalignments(K,y);
talignments,

plot_specs;
plot_imagesc;
%fit_sigmoid;
bar_alignments;
fit_sorted_sigmoid;
fit_diff_sigmoid;
%plot_eigenvalues;

%save_alignments;
%save_plots;

