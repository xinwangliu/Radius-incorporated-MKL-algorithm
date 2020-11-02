tic;
data = load('data/bleakley.mat');
time_data=toc;
struct2workspace(data);   %  K  y  beta_true

alignments = compute_alignments(K);
specs = compute_specs(alignments);
talignments = compute_targetalignments(K,[],Y);
talignments,

%plot_specs;
plot_imagesc;
%fit_sigmoid;
bar_alignments_sorted;
bar_alignments_unsorted;
%fit_sorted_sigmoid;
%fit_diff_sigmoid;
%plot_eigenvalues;

%save_alignments;
%save_plots;

