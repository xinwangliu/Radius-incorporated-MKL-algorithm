function [w1,w2] = mybinarymkl (M, tr_label, pathfile, resultfile, EXPname, permind)
% M: 3 dimensional matrix (n x n x p) of p base kernels on training data
% tr_label: vector of training labels (n x 1)
% pathfile: name of the file that contains various paths for log files, executables, etc
% resultfile: filename for writing the results (text file)
% EXPname: name of the experiment (or dataset) -- is appened to all log files and pegasos files 
% permind (optional): integer that will be appended to all log files and pegasos files (to avoid conflicts among different runs on same dataset) -- is set to '0' if not passed
% w1: kernel weights (obtained when K-classifier is cross-validated on accuracy in K-space
% w2: kernel weights (obtained when K-classifier is cross-validated on hinge-loss in K-space

   if (~exist('permind', 'var'))
      permind = 0;
   end
   
   %%% setup different path variables
   run(pathfile);

   %%% pegasos parameters
   ITER = 100000; BATCH_SIZE = 100;
   SUBSAMPLE = 1;
   lambda_upper = 100;  
   lambda_lower = 0.00000001;  
   divfact = 4; %range of 'lambda' parameter for K-classifier (pegasos) cross-validation, and step size
   SUBSAMPLE = 1;
   subsample_pos = 0;
   subsample_neg = 0;
   valratio = 0.2;  %subsampling factor for validation similarity data
   is_symmetric = 1; %1 if the kernel is symmetric
   %end_pegasos_parameters
   
   num_kernels = size(M,3);
   %%% cross-validation for 'lambda' parameter for K-classifier (pegasos)
   lambda(1) = lambda_upper;
   i = 1;
   while lambda(i)/divfact >= lambda_lower
      lambda(i+1) = lambda(i)/divfact;
      i = i+1;
   end

   MODELFILE = sprintf ('%s/model_%s_perm%d_CV.peg', MODELDIR, EXPname, permind);
   PEGA_LOGFILE = sprintf ('%s/pegasos_%s_perm%d_CV.log', LOGDIR, EXPname, permind);
   sim_fname_fulltr_tr = sprintf ('%s/sim_%s_perm%d_fulltr_tr_valratio%.1f.libsvm', SIMDIR, EXPname, permind, valratio);
   sim_fname_fulltr_ve = sprintf ('%s/sim_%s_perm%d_fulltr_ve_valratio%.1f.libsvm', SIMDIR, EXPname, permind, valratio);
   if (exist(sim_fname_fulltr_tr)~=2 || exist(sim_fname_fulltr_ve)~=2)
      fprintf ('writing files: %s & %s\n', sim_fname_fulltr_tr, sim_fname_fulltr_ve);
      targetmat = construct_targetmat(tr_label, tr_label);
      write_svm_trval(sim_fname_fulltr_tr, sim_fname_fulltr_ve, M, targetmat, is_symmetric, subsample_pos, subsample_neg, SUBSAMPLE, valratio);
   else
      fprintf ('using already existing files %s and %s ...\n', sim_fname_fulltr_tr, sim_fname_fulltr_ve);
   end
   pegasos_command = sprintf ('unset LD_LIBRARY_PATH; %s -lambdaFile %s -k %d -iter %d -modelFile %s -logFile %s -testFile %s %s', PEGASOSMult_EXEC, LAMBDAFILE, BATCH_SIZE, ITER, MODELFILE, PEGA_LOGFILE, sim_fname_fulltr_ve, sim_fname_fulltr_tr)
   system (pegasos_command);
   for i=1:length(lambda)
      pegalog = sprintf ('%s%d', PEGA_LOGFILE, i);
      fid_pegalog = fopen (pegalog, 'r');
      if (fid_pegalog == -1)
	 fprintf ('cant open %s for reading !! \n', pegalog);
      end
      col1 = textscan (fid_pegalog, '%f%*[^\n]');
      simloss_ve(i) = col1{1}(7);
      simacc_ve(i) = 1-col1{1}(8);
      fclose (fid_pegalog);
   end
   %end_pegasos_cross_validation


   %%% run pegasos with selected 'lambda'
   MODELFILE = sprintf ('%s/model_%s_perm%d.peg', MODELDIR, EXPname, permind);
   PEGA_LOGFILE = sprintf ('%s/pegasos_%s_perm%d.log', LOGDIR, EXPname, permind); 
   fid_binmkl = fopen (resultfile, 'w');
   if (fid_binmkl == -1)
      fprintf ('cant open %s for writing !!\n', resultfile);
   end
   fprintf (fid_binmkl, '\n %s (%d) \n\n', EXPname, permind);
   for ii=1:length(lambda)
      fprintf (fid_binmkl, 'simacc: %.5f   simloss: %.5f   lambda: %.9f\n', simacc_ve(ii), simloss_ve(ii), lambda(ii));
   end
   [val ind] = max(simacc_ve);
   lambda_simacc = lambda(ind);
   fprintf (fid_binmkl, '\nmax simacc: %.5f   at   lambda: %.9f\n', val, lambda(ind));
   [val ind] = min(simloss_ve);
   lambda_simloss = lambda(ind);
   fprintf (fid_binmkl, 'min simloss: %.5f   at   lambda: %.9f\n', val, lambda(ind));
   sim_fname = sprintf ('%s/sim_%s_perm%d_fulltr.libsvm', SIMDIR, EXPname, permind);
   if (exist(sim_fname)~=2)
      fprintf ('writing %s for training K-classifier (pegasos) ...\n', sim_fname);
      targetmat = construct_targetmat(tr_label, tr_label);
      write_svm(sim_fname, M, targetmat, is_symmetric, subsample_pos, subsample_neg, SUBSAMPLE);
   else
      fprintf ('using already existing %s to train K-classifier (pegasos)...\n', sim_fname);
   end

   %% run pegasos with 'lambda' selected using cross-validation on accuracy of K-classifier
   pegasos_command = sprintf ('unset LD_LIBRARY_PATH; %s -lambda %.9f -k %d -iter %d -modelFile %s %s > %s', PEGASOS_EXEC, lambda_simacc, BATCH_SIZE, ITER, MODELFILE, sim_fname, PEGA_LOGFILE)
   system (pegasos_command);
   w1 = read_model(MODELFILE, num_kernels);

   %% run pegasos with 'lambda' selected using cross-validation on hinge-loss of K-classifier
   pegasos_command = sprintf ('unset LD_LIBRARY_PATH; %s -lambda %.9f -k %d -iter %d -modelFile %s %s > %s', PEGASOS_EXEC, lambda_simloss, BATCH_SIZE, ITER, MODELFILE, sim_fname, PEGA_LOGFILE)
   system (pegasos_command);
   w2 = read_model(MODELFILE, num_kernels);

   %% writing kernel weights
   fprintf (fid_binmkl, '\n\nkernel weights (K-classifier cross-validated on accuracy in K-space):\n');
   for i=1:length(w1)
      fprintf (fid_binmkl, '%f ', w1(i));
   end
   fprintf (fid_binmkl, '\n\n\nkernel weights (K-classifier cross-validated on hinge loss in K-space):\n');
   for i=1:length(w2)
      fprintf (fid_binmkl, '%f ', w2(i));
   end
