function [trST,trSTp,xsup,w,d,pos,alpha,obj]=traceSTclass(x,y,C,lambda,K0,Sigma,verbose,alphainit)
 
n=size(K0,1);

%compute the radius of enclosing minimum ball
[trST, trSTp]= traceST(K0,Sigma);		

%----------------------------------------------------------------------
%      monqp(H,b,c) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x <= c 
%     x    
%----------------------------------------------------------------------
K= sumKbeta(K0,Sigma);

H=(y*y').*K;
H = H/trST;

f = ones(n,1);
A=y;
b=0;
U=C*ones(n,1);

if verbose ~= 0 
    disp('in QP'); 
end;

[alpha , b0, pos] = monqp(H,f,A,b,U,lambda,verbose,x,K,alphainit);  

if verbose ~= 0 
    disp('out QP'); 
end;

alphaall=zeros(n,1);
alphaall(pos)=alpha;
obj=-0.5*alphaall'*H*alphaall +f'*alphaall;

if ~isempty(x)
    xsup=x(pos,:);
else
    xsup=[];
end
ysup=y(pos);
w=(alpha.*ysup);
d=b0;