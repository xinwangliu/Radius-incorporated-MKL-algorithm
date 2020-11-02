function result = exp_toy_motoaki_marius

%--------------------------------------
%  Parameters
%--------------------------------------

plot_gap = 1;

%gaps = 10.^[-10, -1.35, -0.9, -0.5, -0.2 0.5];
gaps =fliplr( [1, 4, 9, 17, 28, 50]);
ntrs = 50;
nte = 10000;
nval = 10000;
reps = 250;
Cs = logspace(-2.5,-0.5,6);
N = 50;
d =1;
sep = 3.5;
verbosity = 0;
epsilon = 10^(-3);
ridge = 10^(-10);
mkl_norms = [1, 4/3, 2, 4, inf];
duality_gap = 0.05;

%--------------------------------------
% Assertions & Cases
%--------------------------------------
tic;
if plot_gap
  assert(max(size(ntrs))==1);
end
n = max(ntrs) + nval + nte;
assert(mod(n,2)==0);
assert(all(mod(ntrs,2)==0));
assert(all(gaps>0));
% 
%--------------------------------------
% Main Loops
%--------------------------------------

tic; t=0;
for a=1:length(gaps)
  gap = gaps(a);

  %-------------------
  % Prepare mu's
  %-------------------
  block_mu = [ones(1,gap) zeros(1,N-gap)]';
  block_mu = block_mu/max(abs(block_mu));
  mu = fvnorm(reshape(repmat(block_mu',d,1),N,1));
  mu = sep * mu / norm(mu);
  block_mus(:,a) = block_mu;
  %mus(:,a) = mu;

  for r =1:reps
    %-------------------
    % generate toy data
    %-------------------
    X =[ randn(N,max(ntrs)+nval+nte) + repmat( [0.5*mu,-0.5*mu] ,1 ,(max(ntrs)+nval+nte)/2 ) ];
    y= repmat([-1 1],  1, (max(ntrs)+nval+nte)/2);
    
    %-------------------
    % calculate kernels for gap/ntr experiment
    %-------------------
    X=reshape(X', [1 fliplr(size(X))]);
    for ni = 1:length(ntrs)
      %-------------------
      % divide
      %-------------------
      ntr = ntrs(ni);
      divtr = [1:ntr];
      divval = [max(ntrs)+1:max(ntrs)+nval];
      divte = [max(ntrs)+nval+1:max(ntrs)+nval+nte];

      for m = 1:length(mkl_norms)
        mkl_norm = mkl_norms(m);
        for c = 1:length(Cs)
          %-------------------
          % train
          %-------------------
          C = Cs(c);
          t=t+1;
          svm = train_sgmklSVM( 'X', X(1,divtr,:), 'y', y(1,divtr), 'C',C,'verbosity',verbosity, 'mkl_norm',mkl_norm,'epsilon',epsilon,'normalization','VARIANCE','kernel','linear','mkl_solver','DIRECT', 'cachesize', 5, 'duality_gap',duality_gap);
          duality_gaps(a,m,c,r,ni) = svm.duality_gap;
          outval = apply_sgmklSVM( 'svm', svm, 'X', X(1,divval,:),'verbosity',verbosity, 'mkl_norm', mkl_norm,'epsilon',epsilon);
          outte = apply_sgmklSVM( 'svm', svm, 'X', X(1,divte,:),'verbosity',verbosity, 'mkl_norm', mkl_norm, 'epsilon',epsilon');
          errs_val(a,ni,m,c,r) = mean( y(1,divval) ~= sign(outval) );
          errs(a,ni,m,c,r) = mean( y(1,divte) ~= sign(outte) );
          delta_beta(a,ni,m,c,r) = norm( svm.beta'/norm(svm.beta) - block_mu/norm(block_mu) );
          print_progress(t, length(gaps)*length(ntrs)*length(mkl_norms)*length(Cs)*reps);
        end
        
      end
      
    end
    
  end

end

%--------------------------------------
% Collecting Results
%--------------------------------------

result.verbose.err=errs(:,1,:,:,:);
stds_val = std(errs_val,1,5);
stderrs_val = stds_val ./ reps;
errs_val = mean(errs_val,5);
[foo c] = min( errs_val + stderrs_val, [], 4 );

% test errors and stds
stds = std(errs,1,5);
stderrs = stds ./ reps;
errs = mean(errs,5);
for i=1:size(c,1)
  for j=1:size(c,2)
    for k=1:size(c,3)
      err(i,j,k) = errs(i,j,k,c(i,j,k));
      sstd(i,j,k) = stds(i,j,k,c(i,j,k));
      stderr(i,j,k) = stderrs(i,j,k,c(i,j,k));
    end
  end
end
result.err = err;
result.std = sstd;
result.stderr = stderr;

% delta beta
stds = std(delta_beta,1,5);
stderrs = stds ./ reps;
errs = mean(delta_beta,5);
for i=1:size(c,1)
  for j=1:size(c,2)
    for k=1:size(c,3)
      beta_err(i,j,k) = errs(i,j,k,c(i,j,k));
      beta_sstd(i,j,k) = stds(i,j,k,c(i,j,k));
      beta_stderr(i,j,k) = stderrs(i,j,k,c(i,j,k));
    end
  end
end

result.beta_err = beta_err;
result.beta_std = beta_sstd;
result.beta_stderr = beta_stderr;

result.bayes_err = 1-normal_cdf(sep/2);
result.err_indices = ['gaps/', 'ntrs', 'mkl_norms', 'Cs', 'reps'];
result.models.errs = errs;
result.models.stds = stds;
result.models.stderrs = stderrs;
result.C = Cs(c);
result.Cs = Cs;
result.block_mus =  block_mus;

result.paras.gaps = gaps;
result.paras.ntrs = ntrs;
result.paras.mkl_norms = mkl_norms;
result.paras.Cs = Cs;
result.paras.reps = reps;
result.paras.nte = nte;
result.paras.sep = sep;
result.reps = reps;
result.N=N;
result.d=d;
duality_gaps = max(duality_gaps,[],4);
result.duality_gap = max(duality_gaps,[],5);
result.ridge = ridge;

%--------------------------------------
% gap Plot
%--------------------------------------

LineWidth = 1.5;
FontSize = 12;

if plot_gap
  %--------------------------------------
  % true model sparsity
  %--------------------------------------
  figure(1);
  clf;
  l = floor(sqrt(length(gaps)))+1;
  for p=1:length(gaps)
    subplot(l,l,p);
    bar(block_mus(:,p)/max(block_mus(:,p)));
    %axis([0 6 0 1]);
    axis off;
  end
  %--------------------------------------
  % test  errors
  %--------------------------------------
  figure(2);
  clf;
  hold on;
  title(['TOY-Motoaki:    ntr=' sprintf('%0.4g',ntrs) '   nte=' sprintf('%0.4g',nte) , ...
             '    N=' int2str(N)  '     d=' int2str(d)   '     rep=' int2str(reps)  , ...
                '     no.Cs=' int2str(length(Cs))] );
  h=xlabel('data sparsity');
  set(h,'FontSize',FontSize);
  h=ylabel('test error');
  set(h,'FontSize',FontSize);
  h1=errorbar([1:length(gaps)], squeeze(err(:,1,1)), squeeze(sstd(:,1,1)),'r');
  h2=errorbar([1:length(gaps)], squeeze(err(:,1,2)), squeeze(sstd(:,1,2)),'m');
  h3=errorbar([1:length(gaps)], squeeze(err(:,1,3)), squeeze(sstd(:,1,3)),'b');
  h4=errorbar([1:length(gaps)], squeeze(err(:,1,4)), squeeze(sstd(:,1,4)),'g');
  h5=errorbar([1:length(gaps)], squeeze(err(:,1,5)), squeeze(sstd(:,1,5)),'k');
  h6 = plot([1:length(gaps)], result.bayes_err*ones(1,length(gaps)),':k');
  hs = [h1(2) h2(2) h3(2) h4(2) h5(2) ];
  %hs = [h1(2) h2(2) h3(2) ];
   for h = hs
    set(h,'LineWidth',LineWidth);
  end
  set(h6,'LineWidth',LineWidth+1);
  for i=2:length(mkl_norms)-1
    str{i-1} = sprintf('%0.3g',mkl_norms(i));
  end
  %hh =legend(hs ,'L_1',['L_{', str{1}, '}'],['L_{', str{2}, '}'],['L_{', str{3}, '}'], 'L_{\infty}','bayes error');
  hh =legend(hs ,'L_1',['L_{', str{1}, '}'], 'L_{\infty}');
  set(hh,'FontSize',FontSize);
  %--------------------------------------
  % beta errors
  %--------------------------------------
  figure(3);
  clf;
  hold on;
  title(['TOY-MOTOAKI:    ntr=' sprintf('%0.4g',ntrs) '   nte=' sprintf('%0.4g',nte) , ...
             '    N=' int2str(N)  '     d=' int2str(d)   '     rep=' int2str(reps)  , ...
                '     no.Cs=' int2str(length(Cs))] );
  h=xlabel('data sparsity');
  set(h,'FontSize',FontSize);
  h=ylabel('\Delta\theta');
  set(h,'FontSize',FontSize);
  h1=errorbar([1:length(gaps)], squeeze(beta_err(:,1,1)), squeeze(beta_sstd(:,1,1)),'r');
  h2=errorbar([1:length(gaps)], squeeze(beta_err(:,1,2)), squeeze(beta_sstd(:,1,2)),'m');
  h3=errorbar([1:length(gaps)], squeeze(beta_err(:,1,3)), squeeze(beta_sstd(:,1,3)),'b');
  h4=errorbar([1:length(gaps)], squeeze(beta_err(:,1,4)), squeeze(beta_sstd(:,1,4)),'g');
  h5=errorbar([1:length(gaps)], squeeze(beta_err(:,1,5)), squeeze(beta_sstd(:,1,5)),'k');
  hs = [h1(2) h2(2) h3(2) h4(2) h5(2)];
  %hs = [h1(2) h2(2)  h3(2) ];
   for h = hs
    set(h,'LineWidth',LineWidth);
  end
  for i=2:length(mkl_norms)-1
    str{i-1} = sprintf('%0.3g',mkl_norms(i));
  end
  %hh =legend(hs ,'L_1',['L_{', str{1}, '}'],['L_{', str{2}, '}'],['L_{', str{3}, '}'], 'L_{\infty}','bayes error');
  hh =legend(hs ,'L_1',['L_{', str{1}, '}'], 'L_{\infty}');
  set(hh,'FontSize',FontSize);
  %--------------------------------------
  %  computed  model sparsity
  %--------------------------------------
   %%%%%%%% NOT IMPLEMENTED
 %save('~/temp/toy_result_gap.mat', 'result');
end

result.time =toc;
result.mkl_norm = mkl_norms;
if length(ntrs)==1
  result.n = [ntrs, nval, nte];
else
  result.ntrs=ntrs;
  result.nvt = [nval nte];
end
result.paras.max_duality_gap = duality_gap;