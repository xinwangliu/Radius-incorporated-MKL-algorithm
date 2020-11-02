function  radiusp = constraintedBasicKernelRadius(K,option)

n=size(K,1);
d= size(K,3);
lambda0 = option.lambdareg;
verbose = option.verbosesvm;

radiusp = zeros(d,1);
for p = 1:d
    H= 2 * K(:,:,p);
    f= diag(K(:,:,p));
    A=ones(n,1);
    b=1;
    B=ones(n,1);
    
    [beta, lambdat, pos] = monqp(H,f,A,b,B,lambda0,verbose);
    
    optimalBeta=zeros(n,1);
    optimalBeta(pos)=beta;
    
    radiusp(p) = -0.5*optimalBeta'*H*optimalBeta + f'*optimalBeta;
end 
radiusp=radiusp(:)';