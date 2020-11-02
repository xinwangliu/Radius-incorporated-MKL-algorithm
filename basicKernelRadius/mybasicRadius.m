function [radiusp]= mybasicRadius(K)

d=size(K,3);
n=size(K,1);

radiusp = zeros(d,1);
%radius = 0;
for p = 1:d
    H= 2 * K(:,:,p);
    f= diag(K(:,:,p));
    A=ones(n,1);
    b=1;
    B=ones(n,1);
    lambda0=1e-10;
    verbose=0;
    [beta, lambdat, pos] = monqp(H,f,A,b,B,lambda0,verbose);
    
    optimalBeta=zeros(n,1);
    optimalBeta(pos)=beta;
    
    radiusp(p) = -0.5*optimalBeta'*H*optimalBeta + f'*optimalBeta;
    %radius = radius +Sigma(p)*radiusp(p);
end 
radiusp=radiusp(:);
