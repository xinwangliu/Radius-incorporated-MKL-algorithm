function specs = compute_specs(alignments);

aux = eig(alignments);
specs = flipud(aux(:));
