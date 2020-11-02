n=500;
div=randperm(2000);
div=div(1:n);

tic;
data = load('data/mklbio.mat');
time_data=toc;
struct2workspace(data);   %  K  y  beta_true
K = K(div,div,:);
y = y(div);
K = kscale(K);


alignments = compute_alignments(K);
specs = compute_specs(alignments);
talignments = compute_targetalignments(K,y);
talignments,

%plot_specs;
plot_imagesc;
%fit_sigmoid;
bar_alignments_unsorted;
bar_alignments_sorted;
%fit_sorted_sigmoid;
%fit_diff_sigmoid;
%plot_eigenvalues;

%save_alignments;
save_plots;
