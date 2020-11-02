function [radius, optimalBeta]= l2SVMradius(K,C,lambda,verbose)

n=size(K,1);
H=2*( K + eye(size(K))/C);
f=(diag(K)+1/C);
A=ones(n,1);
b=1;
B=ones(n,1);

[beta, lambda0, pos] = monqp(H,f,A,b,B,lambda,verbose);

optimalBeta=zeros(n,1);
optimalBeta(pos)=beta;

radius = -0.5*optimalBeta'*H*optimalBeta + f'*optimalBeta;

% 
