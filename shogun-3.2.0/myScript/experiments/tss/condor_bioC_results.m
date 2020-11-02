function res = condor_bioC_results

fprintf('for full Plots with logistic errorbars outcomment last lines of m file\n');

% hier die Parameter aus condor_bioC herkopieren
ntrain = [ 1250 2500 5000 10000 20000 40000 60000 ]; %max 60000
ntrain=1250;
nr_test = 4;  % the index number of the test set partition  (1-4)
Cs = [ 0.25  1 2 4 8 16 32];
mkl_norms = [0, 1.333, 1,  2, 4];
mkl_norms =2


dirname = ['~/condor/condor_bioC'];
cd(dirname);

%-----------collect results for normed attack
nonsparse_norms = setdiff(mkl_norms, [0,1]);
ntr=ntrain;
n=ntr;

for mkl_norm=mkl_norms
  for  i=1:length(ntr)
    n=ntr(i);
    reps(i) = floor(60000/n);
    for j=1:length(Cs)
      for rep = [1:reps(i)]

        C=Cs(j);

        %err_filename = ['n=', int2str(n), '/C=', num2str(C), ...
        %  '/mkl_norm=',int2str_leading_zeros(mkl_norm,1), '/rep=', int2str(rep), '/stderr'];
        %assert(exist(err_filename, 'file'));
        %e = load(err_filename);
        %if ~max(max(size(e)))  == 0
        %  e,
        %  warning('Error during in stderr ', err_filename);
        %end
        filename = ['nr_test=', int2str(nr_test),'ntr=', int2str(n), 'mkl_norm=', sprintf('%0.4g',mkl_norm), ...
                          'C=',  sprintf('%0.4g',C), 'rep=', int2str_leading_zeros(rep,2), '.mat'];
        if ~(exist(filename, 'file')),
          fprintf('\nFile %s does not exist.\n', filename);
          assert(exist(filename, 'file'));
        end
        fprintf('%d ',rep);
        load(filename);
        l=load(filename);
        duality_gap = condor_bioC_check_duality_gap(l.divtr,l.alphas,l.b,l.w,mkl_norm);
          keyboard;
        if mkl_norm==0
          av.ntr = ntr;
          avROC_val{i,j}(rep) = auROC_val;
          avPRC_val{i,j}(rep) = auPRC_val;
          avROC_te{i,j}(rep) = auROC_te;
          avPRC_te{i,j}(rep) = auPRC_te;
          avBeta{i,j}(rep,:) = l.w;
        elseif mkl_norm==1
          l1.ntr = ntr;
          l1ROC_val{i,j}(rep) = auROC_val;
          l1PRC_val{i,j}(rep) = auPRC_val;
          l1ROC_te{i,j}(rep) = auROC_te;
          l1PRC_te{i,j}(rep) = auPRC_te;
          l1Beta{i,j}(rep,:) = l.w;
        elseif mkl_norm>1
          k = find(mkl_norm==nonsparse_norms);
          l2{k}.ntr = ntr;
          l2ROC_val{i,j,k}(rep) = auROC_val;
          l2PRC_val{i,j,k}(rep) = auPRC_val;
          l2ROC_te{i,j,k}(rep) = auROC_te;
          l2PRC_te{i,j,k}(rep) = auPRC_te;
          l2Beta{i,j,k}(rep,:) = l.w;
      end
    end
  end
 end
end

for  i=1:length(ntr)
  for  j=1:length(Cs)
    av.val.auROC_std(i,j) = std(avROC_val{i,j},1);
    av.val.auROC_stderr(i,j) = av.val.auROC_std(i,j) / sqrt(size(avROC_val{i,j},2));
    av.val.auROC(i,j) = mean(avROC_val{i,j});
    l1.val.auROC_std(i,j) = std(l1ROC_val{i,j},1);
    l1.val.auROC_stderr(i,j) = l1.val.auROC_std(i,j) / sqrt(size(l1ROC_val{i,j},2));
    l1.val.auROC(i,j) = mean(l1ROC_val{i,j});
    for k=1:length(nonsparse_norms)
      l2{k}.mkl_norm = nonsparse_norms(k);
      l2{k}.val.auROC_std(i,j) = std(l2ROC_val{i,j,k},1);
      l2{k}.val.auROC_stderr(i,j) = l2{k}.val.auROC_std(i,j) / sqrt(size(l2ROC_val{i,j,k},2));
      l2{k}.val.auROC(i,j) = mean(l2ROC_val{i,j,k});
  %--------
%    av.val.auPRC_std(i,:) = std(avPRC_val{i},[],1);
%    av.val.auPRC_stderr(i,:) = av.val.auPRC_std(i,:) / sqrt(size(avPRC_val{i},2));
%    av.val.auPRC(i,:) = mean(avPRC_val{i},1);
%    [foo, ind] = max(mean(avPRC_val{i},1)-av.val.auPRC_stderr(i,:));
%    av.C_prc(i) = Cs(ind);
%    av.auPRC(i) = mean(avPRC_te{i}(:,ind),1);
%    av.auPRC_std(i) = std(avPRC_te{i}(:,ind),[],1);
%    av.auPRC_stderr(i) = std(avPRC_te{i}(:,ind),[],1) / sqrt(length(avPRC_te{i}(:,ind))) ;
%    l1.val.auPRC_std(i,:) = std(l1PRC_val{i},[],1);
%    l1.val.auPRC_stderr(i,:) = l1.val.auPRC_std(i,:) / sqrt(size(l1PRC_val{i},2));
%    l1.val.auPRC(i,:) = mean(l1PRC_val{i},1);
%    [foo, ind] = max(mean(l1PRC_val{i},1)-av.auPRC_stderr(i));
%    l1.C_prc(i) = Cs(ind);
%    l1.auPRC(i) = mean(l1PRC_te{i}(:,ind),1);
%    l1.auPRC_std(i) = std(l1PRC_te{i}(:,ind),[],1);
%    l1.auPRC_stderr(i) = std(l1PRC_te{i}(:,ind),[],1) / sqrt(length(l1PRC_te{i}(:,ind))) ;
%    l2{k}.val.auPRC_std(i,:) = std(l2PRC_val{i,k},[],1);
%    l2{k}.val.auPRC_stderr(i,:) = l2{k}.val.auPRC_std(i,:) / sqrt(size(l2PRC_val{i,k},2));
%    l2{k}.val.auPRC(i,:) = mean(l2PRC_val{i,k},1);
%    [foo, ind] = max(mean(l2PRC_val{i,k},1)-l2{k}.val.auPRC_stderr(i,:));
%    l2{k}.C_prc(i) = Cs(ind);
%    l2{k}.auPRC(i) = mean(l2PRC_te{i,k}(:,ind),1);
%    l2{k}.auPRC_std(i) = std(l2PRC_te{i,k}(:,ind),[],1);
%    l2{k}.auPRC_stderr(i) = std(l2PRC_te{i,k}(:,ind),[],1) / sqrt(length(l2PRC_te{i,k}(:,ind))) ;
    end
  end
  [foo, j] = max( av.val.auROC(i,:) -  av.val.auROC_stderr(i,:) );
  av.C_roc(i) = Cs(j);
  av.auROC(i) = mean(avROC_te{i,j});
  av.auROC_std(i) = std(avROC_te{i,j},1);
  av.auROC_stderr(i) = av.auROC_std(i) / sqrt(length(avROC_te{i,j})) ;
  av.beta(i,:)= mean(avBeta{i,j},1)
  [foo, j] = max(l1.val.auROC(i,:) -  l1.val.auROC_stderr(i,:));
  l1.C_roc(i) = Cs(j);
  l1.auROC(i) = mean(l1ROC_te{i,j});
  l1.auROC_std(i) = std(l1ROC_te{i,j},1);
  l1.auROC_stderr(i) = l1.auROC_std(i) / sqrt(length(l1ROC_te{i,j})) ;
  l1.beta(i,:)= mean(l1Beta{i,j},1)
  [foo, j] = max( av.val.auROC(i,:) -  av.val.auROC_stderr(i,:) );
  for k=1:length(nonsparse_norms)
    l2{k}.C_roc(i) = Cs(j);
    l2{k}.auROC(i) = mean(l2ROC_te{i,j,k});
    l2{k}.auROC_std(i) = std(l2ROC_te{i,j,k},1);
    l2{k}.auROC_stderr(i) = l2{k}.auROC_std(i) / sqrt(length(l2ROC_te{i,j,k})) ;
    l2{k}.beta(i,:)= mean(l2Beta{i,j,k},1)
  end
end

% collect betas
beta(1,:)=l1.beta(end,:);
for k=1:length(nonsparse_norms)
  beta(k+1,:)=l2{k}.beta(end,:);
end
title(int2str(ntr(end)));
beta(end+1,:)=av.beta(end,:);

%plot betas
figure;
for i=1:5
  subplot(2,3,i);
  bar(beta(i,:));
end

%  figure;
%  hold on;
%  title('training set size vs. area under PR curve')
%  h1=errorbar(ntr, av.auPRC, av.auPRC_stderr,'k');
%  h2=errorbar(ntr, l1.auPRC, l1.auPRC_stderr,'r');
%  h3=errorbar(ntr, l2{1}.auPRC, l2{1}.auPRC_stderr,'b-');
%  h4=errorbar(ntr, l2{2}.auPRC, l2{2}.auPRC_stderr,'b:');
%  handle = [h1(2) h2(2) h3(2) h4(2)];
%  legend(handle ,'canonical','l_1',['l_' str{1}], ['l_' str{2}]);

for k=1:length(nonsparse_norms)
  str{k} = sprintf('%0.4g',nonsparse_norms(k));
end

figure;
hold on;
title('training set size vs. area under ROC curve');
plot(ntr(end-2:end), l1.auROC(end-2:end),'r');
plot(ntr(end-2:end), l2{1}.auROC(end-2:end),'m');
plot(ntr(end-2:end), l2{2}.auROC(end-2:end), 'b');
plot(ntr(end-2:end), l2{3}.auROC(end-2:end), 'g');
plot(ntr(end-2:end), av.auROC(end-2:end),'k');
h1=errorbar(ntr(1:end-2), l1.auROC(1:end-2), l1.auROC_stderr(1:end-2),'r');
h2=errorbar(ntr(1:end-2), l2{1}.auROC(1:end-2), l2{1}.auROC_stderr(1:end-2),'m');
h3=errorbar(ntr(1:end-2), l2{2}.auROC(1:end-2), l2{2}.auROC_stderr(1:end-2),'b');
h4=errorbar(ntr(1:end-2), l2{3}.auROC(1:end-2), l2{3}.auROC_stderr(1:end-2),'g');
h5=errorbar(ntr(1:end-2), av.auROC(1:end-2), av.auROC_stderr(1:end-2),'k');
handle = [h1(2) h2(2) h3(2) h4(2)];
legend(handle ,'L_1',['L_{', str{1}, '}'], ['L_{', str{2}, '}'],  ['L_{', str{3}, '}'], 'L_{\infty}');

res.av=av;
res.l1=l1;
res.l2=l2;
res.ntrain = ntrain;
res.Cs = Cs;
res.nr_test = nr_test;
res.mkl_nornms = mkl_norms;

%cd ~;
%save bioC;
% FÜR VOLLSTÄNDIGE PLOTS WEGEN ERRORBAR LOG PROBLEM
% DATEI "bioC" nach E:\ kopieren und bioC_results AUF NOTEBOOK AUFRUFEN

