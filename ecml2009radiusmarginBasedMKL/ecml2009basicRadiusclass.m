function [w,d,pos,obj]= ecml2009basicRadiusclass(K0,Sigma,y,C,radiusp,option)
 
n=size(K0,1);

verbose = option.verbosesvm;
lambda0 = option.lambdareg;
%compute the radius of enclosing minimum ball		
radius  = radiusp' * Sigma;
%----------------------------------------------------------------------
%      monqp(H,b,c) solves the quadratic programming problem:
% 
%    min 0.5*x'Hx - d'x   subject to:  A'x = b  and  0 <= x <= c 
%     x    
%----------------------------------------------------------------------
K = sumKbeta(K0,Sigma);
H = (y*y').*K + eye(n)* radius/C;

f = ones(n,1);
A=y;
b=0;

[alpha , b0, pos] = monqpCinfty(H,f,A,b,lambda0,verbose);  

alphaall=zeros(n,1);
alphaall(pos)=alpha;
obj=-0.5*alphaall'*H*alphaall +f'*alphaall;

ysup=y(pos);
w=(alpha.*ysup);
d=b0;
