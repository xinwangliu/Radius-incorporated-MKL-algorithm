clear
clc
warning off;

path = 'F:\work2015\';
addpath(genpath(path));
dataName = 'heart'; %%% flower17; flower102; CCV; caltech101_numofbasekernel_10
%% %% washington; wisconsin; texas; cornell
%% caltech101_nTrain5_48; proteinFold
load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
% load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numclass = length(unique(Y));
numker = size(KH,3);
num = size(KH,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KH = kcenter(KH);
KH = knorm(KH);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qnorm = 2;
C = 100;
verbose = 1;
% [H_normalized,gamma,obj] = mkkmeans_train(KH,numclass,qnorm);
% res_gnd = myNMIACC(H_normalized,Y,numclass);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsionset = [0.5:0.1:0.5];
for ie = 1:length(epsionset)
     for iter = 1:1
        load([path,'work2016\generateAbsentMatrix\',dataName,'_missingRatio_',num2str(epsionset(ie)),...
            '_missingIndex_iter_',num2str(iter),'.mat'],'S');
        
% %         %%%%%%%%%%%--Zero-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
        KH1 = algorithm2(KH,S);
        [Sigma1,Alpsup1,w01,pos1,obj1] = mklGLasso_alt0_class(KH1,Y,C,qnorm,verbose);
        Kmatrix1  = mycombFun(KH1,Sigma1);
        ypred1 = mysvmval(Alpsup1,w01,pos1,Kmatrix1);
        acc(1) = mean(sign(ypred1)==Y);

        %%%%%%%%%%%--mean-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
        KH2 = algorithm3(KH,S);
        [Sigma2,Alpsup2,w02,pos2,obj2] = mklGLasso_alt0_class(KH2,Y,C,qnorm,verbose);
        Kmatrix2  = mycombFun(KH2,Sigma2);
        ypred2 = mysvmval(Alpsup2,w02,pos2,Kmatrix2);
        acc(2) = mean(sign(ypred2)==Y);
        
%         %%%%%%%%%%--knn-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         KH3 = algorithm0(KH,S,7);
%         [Sigma3,Alpsup3,w03,pos3,obj3] = mklGLasso_alt0_class(KH3,Y,C,qnorm,verbose);
%         Kmatrix3  = mycombFun(KH3,Sigma3);
%         ypred3 = mysvmval(Alpsup3,w03,pos3,Kmatrix3);
%         acc(3) = mean(ypred3==Y);
        
        [KH4,Sigma4,Alpsup4,w04,pos4,obj4] = myabsentmultikernelbinaryclassification(KH,Y,C,S,qnorm);
        Kmatrix4  = mycombFun(KH4,Sigma4);
        ypred4 = mysvmval(Alpsup4,w04,pos4,Kmatrix4);
        acc(4) = mean(sign(ypred4)==Y);
       end
end