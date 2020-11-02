function globals = exp_setup

% random number generator
globals.rnd_seed = 12345;

% accuracy
globals.eps  = 10^-3;

% svm trade-off
globals.svmC = 1;

% nof repetitions:
globals.repetitions = 5;

% nof kernels/sample size
globals.nofks = [50];
globals.nofex = [100,250,500,1000,2000,4000,8000];
% additional sample sizes for larger scales
globals.largescale = [100,250,500,1000,2000,4000,8000,10000,20000,40000,60000];

% SVM cache
globals.cachesize = 50;

% kernel basis
globals.basis = 1.2;

% ridge
globals.ridge_eps = 0;

return;
