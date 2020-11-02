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
globals.nofex = [1000];     
globals.nofks = [10, 20, 50, 100, 200, 500, 1000];      
% additional kernel sizes for larger scales
globals.largescale = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]; 

% SVM cache
globals.cachesize = 50;

% kernel basis
globals.basis = 1.2;

% ridge
globals.ridge_eps = 0;

return;
