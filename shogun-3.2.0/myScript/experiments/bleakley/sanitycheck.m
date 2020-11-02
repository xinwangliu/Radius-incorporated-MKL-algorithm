
function sanitycheck

disp('Bleakley et al. / sanity check');
disp('Reproducing Table 1 (local)');
disp('------------------------------');

% load adjacency matrix
load adj_metabolic.mat;
A = adj_metabolic;
A(A==0) = -1;

%------------------------------
% EXP
fprintf('EXP\n');
load K_metabolic_exp.mat;
K = K_metabolic_exp;
clear K_metabolic_exp;
auc_exp = svm(K,A);

%-------------------------------
% LOC
fprintf('LOC\n');
load K_metabolic_loc.mat;
K = K_metabolic_loc;
clear K_metabolic_loc;
auc_loc = svm(K,A);


%-------------------------------
% PHY
load K_metabolic_phy.mat;
fprintf('PHY\n');
K = K_metabolic_phy;
clear K_metabolic_phy;
auc_phy = svm(K,A);

%-------------------------------
% INT
load K_metabolic_int.mat;
fprintf('INT\n');
K = K_metabolic_int;
clear K_metabolic_int;
auc_int = svm(K,A);

save 'table1.mat' auc_exp auc_loc auc_phy auc_int;

return;


