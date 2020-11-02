function [w,b,nbsv,pos,obj] = mySVMmulticlassoneagainstall(y,C,K,Sigma)

nbclass = length(unique(y));
% Xinwang Liu
KC = sumKbeta(K,Sigma);
% 3D matrices can not be used as numebr of SV changes
w=[];
b=[];
pos=[];
nbsv=zeros(1,nbclass);
obj=0;
for i=1:nbclass
    yone=(y==i)+(y~=i)*-1;
    %%% [w,d,pos,obj,alpha,timeps]=mySVMclass(y,c,K,lambda,verbose)
    [waux,baux,posaux,objaux]=mySVMclass(yone,C,KC);
    nbsv(i)=length(posaux);
    w   = [w;waux];
    b   = [b;baux];
    pos = [pos;posaux];
    obj = obj+objaux;
end;
