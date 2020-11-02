function [Sigma, Alpsup,w0,pos,nbsv,obj] = mklGLasso_multiclass(K,y,C,qnorm,verbose)

nbkernel=size(K,3);
Sigma = ones(nbkernel,1)*(1/nbkernel)^(1/qnorm);
%------------------------------------------------------------------------------%
% Initialize
%------------------------------------------------------------------------------%
nloop = 1;
loop = 1;
nloopmax = 100;
% initial a solution for future use
Kmatrix = sumKbeta(K,Sigma);
%------------------------------------------------------------------------------%
% Update Main loop
%------------------------------------------------------------------------------%
while loop
    [Alpsup,w0,nbsv,pos,obj(nloop)] = mySVMmulticlassoneagainstall(y,C,Kmatrix);
    %-----------------------------------------
    % Update weigths Sigma
    %----------------------------------------
    Sigmaold = Sigma;
    % note!!!
    Sigma = updateMySigma(Alpsup,nbsv,pos,K,Sigmaold,qnorm);
    %%----------------------------------------------------
    %% Run SVM
    %%----------------------------------------------------
    Kmatrix = sumKbeta(K,Sigma);
%     [Alpsup,w0,nbsv,pos,obj(nloop)] = mySVMmulticlassoneagainstall(y,C,Kmatrix);
%     
    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  (max(abs(Sigma - Sigmaold))<1e-4 || ( nloop>2 && (abs((obj(nloop-1)-obj(nloop))/obj(nloop))...
            <1e-4)|| nloop >nloopmax))
        loop = 0;
        fprintf(1,'variation convergence criteria reached \n');
    end
    nloop = nloop+1; 
end