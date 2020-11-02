%plots imagesc of alignments
figure(2);
clf;

imagesc(alignments-diag(ones(1,size(K,3))));
imagesc(alignments);
colorbar;
hold on;
%brpt = data.no_inf(i);
%ma = max(axis);
%mi = min(axis);
%plot( [brpt brpt], [mi ma], 'r-', 'linewidth',2);
%plot( [mi ma], [brpt brpt] , 'r-', 'linewidth',2);

title('pairwise kernel alignments(=dot-products)');
xlabel('kernel-id');
ylabel('kernel-id');
%xlabel('WD / TSS / exon / ang / ener');
%ylabel('WD / TSS / exon / ang / ener');