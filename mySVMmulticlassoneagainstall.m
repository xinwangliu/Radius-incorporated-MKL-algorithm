function [w,b,nbsv,pos,obj]= mySVMmulticlassoneagainstall(y,C,K)

nbclass = length(unique(y));
w=[];
b=[];
pos =[];
nbsv = zeros(1,nbclass);
obj=0;
for i=1:nbclass
    yone = (y==i)+(y~=i)*-1;
    [waux,baux,posaux,objaux] = mySVMclass(yone,C,K);
    nbsv(i)=length(posaux);
    w = [w; waux];
    b = [b; baux];
    pos = [pos; posaux];
    obj = obj + objaux;
end