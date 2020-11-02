basis = 2;
kseq = basis.^[(1:nofks)-1];

%disp(kseq);

idx = randperm(size(Xtr,1));

X = Xtr(idx(1:nofex),:);

ytr(idx(1:nofex))';

K=[];

for i=1:nofks

    sigma = kseq(i);
    K(:,:,i) = rbf(X,sigma);
   
end;
    
