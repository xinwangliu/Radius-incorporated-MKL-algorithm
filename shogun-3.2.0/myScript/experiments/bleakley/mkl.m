
function mkl

% Bleakley exp. -- MKL, 5 fold CV

%----------------------------------
% some params
svm_eps=1e-3;	 % svm epsilon
mkl_eps=1e-3;	 % mkl epsilon
mkl_norm = 2;    % mkl norm
cachesize = 50;
r = 1e-8;

%------------------------------
% load adjacency matrix
load adj_metabolic.mat;
A = adj_metabolic;
A(A==0) = -1;
n = size(A,1);

%------------------------------
% load and perpare kernels
load K_metabolic_exp.mat;
K = K_metabolic_exp + r * trace(K_metabolic_exp) * eye(n);
clear K_metabolic_exp;

load K_metabolic_loc.mat;
K(:,:,2) = K_metabolic_loc + r * trace(K_metabolic_loc) * eye(n);
clear K_metabolic_loc;

load K_metabolic_phy.mat;
K(:,:,3) = K_metabolic_phy + r * trace(K_metabolic_phy) * eye(n);
clear K_metabolic_phy;

C = [1]; % ,10,100,1000];

% init random numbers
rand('seed',26345);
randn('seed',26345);

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
    for i = 1:R

        % prepare matrices
        Ktr = K(folds~=i,folds~=i,:);
        Kte = K(folds~=i,folds==i,:);
    
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
            Kj(j,:,:) = [];
            Kj(:,j,:) = [];
        
            % init MKL
            sg('clean_kernel');
            sg('clean_features', 'TRAIN');
            
            sg('new_classifier', 'MKL_CLASSIFICATION'); 
            sg('mkl_use_interleaved_optimization', 1); % 1 interleaved, 0 wrapper
            %sg('set_constraint_generator', 'LIBSVM');
            sg('set_solver', 'DIRECT'); % DIRECT, NEWTON, CPLEX, AUTO, GLPK
            sg('mkl_parameters', mkl_eps, 0, mkl_norm);
            sg('svm_epsilon', svm_eps);
           
            sg('svm_use_bias',1);
            %sg('c',C(c)*n_pos,C(c)*n_neg);
            sg('c',C(c));
            sg('set_labels','TRAIN',yj');
            sg('set_kernel','COMBINED',cachesize);
            for kk = 1:size(K,3)
                sg('add_kernel',1,'CUSTOM',Kj(:,:,kk),'FULL');
            end;
                      
            % train classifier
            sg('train_classifier');
            
            % classify test set
            sg('clean_kernel');
            sg('clean_features', 'TEST');
            sg('set_labels','TEST',Ate(j,:));
            Kj = Kte;
            Kj(j,:,:) = [];
            sg('set_kernel','COMBINED',cachesize);
            for kk = 1:size(K,3)
                sg('add_kernel',1,'CUSTOM',Kj(:,:,kk),'FULL');
            end;
            Pred(j,:) = sg('classify');
       
        end;
        % evaluate 
        [auc(i)] = rocscore(Pred(:),Ate(:));
    end;
    
    fprintf('%f-norm, C = %f: AUC = %f +- %f\n',mkl_norm,C(c),mean(auc)*100,std(auc)*100/sqrt(R));
    AUC(c,1:3) = [C(c),mean(auc)*100,std(auc)*100/sqrt(R)];
    
end;

str = sprintf('mkl_%f',mkl_norm);
str = strrep(str,'.','_');
str = [str,'.mat'];
save(str,'AUC');

exit;

return;