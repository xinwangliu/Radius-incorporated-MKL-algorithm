clear
clc
warning off;

path = 'F:\work2015\';
addpath(genpath(path));
dataName = 'proteinFold'; %%% flower17; flower102; CCV; caltech101_numofbasekernel_10
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
        [Sigma1,Alpsup1,w01,pos1,nbsv1,obj1] = mklGLasso_multiclass(KH1,Y,C,qnorm,verbose);
        [ypred1,label1,maxi1]= mymultival(Alpsup1,w01,pos1,nbsv1,KH1,Sigma1);
        acc(1) = mean(label1==Y);

        %%%%%%%%%%%--mean-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
        KH2 = algorithm3(KH,S);
        [Sigma2,Alpsup2,w02,pos2,nbsv2,obj2] = mklGLasso_multiclass(KH2,Y,C,qnorm,verbose);
        [ypred2,label2,maxi2]= mymultival(Alpsup2,w02,pos2,nbsv2,KH2,Sigma2);
        acc(2) = mean(label2==Y);
        
%         %%%%%%%%%%--knn-Filling--%%%%%%%%%%%%%%%%%%%%%%%%%%
%         KH3 = algorithm0(KH,S,7);
%         [Sigma3,Alpsup3,w03,pos3,nbsv3,obj3] = mklGLasso_multiclass(KH3,Y,C,qnorm,verbose);
%         [ypred3,label3,maxi3]= mymultival(Alpsup3,w03,pos3,nbsv3,KH3,Sigma3);
%         acc(3) = mean(label3==Y);
        
        [KH4,Sigma4,Alpsup4,w04,pos4,nbsv4,obj4] = myabsentmultikernelclassification(KH,Y,C,S,qnorm);
        [ypred4,label4,maxi4]= mymultival(Alpsup4,w04,pos4,nbsv4,KH4,Sigma4);
        acc(4) = mean(label4==Y);
       end
end