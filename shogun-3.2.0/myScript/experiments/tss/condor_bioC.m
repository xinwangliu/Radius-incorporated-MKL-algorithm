ntrain = [ 1250 2500 5000 10000 20000 40000 60000  ]; %max 20000
nr_test = 4;  % the index number of the test set partition  (1-4)
Cs = [ 0.25  1 2 4 8 16 32 ];
mkl_norms = [0, 1, 1.333, 2, 4];


dir    = '/home/neuro/klomar/condor';      % Where to store files
base   = 'condor_bioC';
script = 'condor_bioC_script';           % The script to execute

addpath /home/neuro/klomar/condor
c = condor(10^6);
size = 500;                              % The memory needed for that

div = load('~/condor/divisions_bioC/divisions.mat');
arguments.divval= div.divval;
div=div.div;
arguments.nr_test = nr_test;

dirname = ['~/condor/condor_bioC'];
cd(dirname);

for nr_te = nr_test
  arguments.divte = div{nr_test};
  nr_train = setdiff( [1:4], nr_test );
  for ntr = ntrain
    for mkl_norm = mkl_norms
      for C=Cs
        dl = length(nr_train);
        reps = floor(20000/ntr);
        if ntr>20000
          dl = 1;
          reps = 1;
          arguments.divtr = [div{1} div{2}];
          if ntr==60000
            arguments.divtr = [div{1} div{2} div{3}];
          end
        end
        for d=1:dl
          nr_tr = nr_train(d);
          for k = 1:reps
            if ntr<=20000
              arguments.divtr = div{nr_tr}((k-1)*ntr+1:k*ntr);
            end
            rep = (d-1)*reps + k;
            arguments.ntr = ntr;
            arguments.rep = rep;
            arguments.C = C;
            arguments.mkl_norm = mkl_norm;
            filename = ['nr_test=', int2str(nr_test),'ntr=', int2str(ntr), 'mkl_norm=', sprintf('%0.4g',mkl_norm), ...
                      'C=',  sprintf('%0.4g',C), 'rep=', int2str_leading_zeros(rep,2), '.mat'];
            if exist(filename, 'file')
              printf('Skipping enqueuing\n');
            else
              job = condor_job(dir, base, script, arguments, size, 'nr_test','ntr','C', 'mkl_norm','rep');
              c = schedule(c, set(job,'priority','6'), 'ignore-all');
            end
          end
        end
      end
    end
  end
end
