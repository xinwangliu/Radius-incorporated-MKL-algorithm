
% Bleakley exp. -- SVM w/ single kernel and 5 fold CV
function AUC = svm(K,A)

C = [1]; % HACIK:::: ,10,100,1000];

% init random numbers
rand('seed',2345);
randn('seed',2345);

% permutation
n = size(K,1);
idx = randperm(n);

% create 5 folds
for i = 1:n
    folds(i) = mod(i,5);
end;
folds(folds==0) = 5;

% nof folds
R = 5;

for c = 1:length(C)

    % 5-fold cv
    for i = 1%:R HACK!!!!

        % prepare matrices
        Ktr = K(folds~=i,folds~=i);
        Kte = K(folds~=i,folds==i);
    
        Atr  = A(folds~=i,folds~=i);
        Ate  = A(folds~=i,folds==i);
        Pred = zeros(size(Ate));
    
        % train local models
        for j = 1:size(Ktr,1)

            % j-th label vector
            yj = Atr(:,j);
            yj(j) = [];
            % j-th kernel matix
            Kj = Ktr;
            Kj(j,:) = [];
            Kj(:,j) = [];
        
            %n_pos = sum(yj == -1)/length(yj);
            %n_neg = sum(yj == +1)/length(yj);
        
            % init svm
            sg('new_svm','LIGHT');
            sg('svm_use_bias',1);
            %sg('c',C(c)*n_pos,C(c)*n_neg);
            sg('c',C(c));
            sg('set_labels','TRAIN',yj');
            sg('set_kernel','CUSTOM',Kj,'FULL');
        
            % train SVM
            sg('train_classifier');
            
            % classify test set
            sg('set_labels','TEST',Ate(j,:));
            Kj = Kte;
            Kj(j,:) = [];
            sg('set_kernel','CUSTOM',Kj,'FULL');
        
            Pred(j,:) = sg('classify');
       
        end;
        % evaluate 
        [auc(i)] = rocscore(Pred(:),Ate(:));
    end;
    
    fprintf('C = %f: AUC = %f +- %f\n',C(c),mean(auc)*100,std(auc)*100/sqrt(R));
    AUC(c,1:3) = [C(c),mean(auc)*100,std(auc)*100/sqrt(R)];
    
end;


