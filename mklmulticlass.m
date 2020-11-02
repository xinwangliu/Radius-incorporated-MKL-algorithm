function [Sigma,Alpsup,w0,pos,nbsv,SigmaH,obj] = mklmulticlass(K,yapp,C,option,verbose)

nbkernel=size(K,3);
Sigma=ones(1,nbkernel)/nbkernel;
loop = 1;
%------------------------------------------------------------------------------%
% Initialize
%------------------------------------------------------------------------------%
Kmatrix=sumKbeta(K,Sigma);

[Alpsup,w0,nbsv,pos,obj0]=mySVMmulticlassoneagainstall(yapp,C,Kmatrix);
[grad] = myGradsvmoneagainstall(K,pos,Alpsup,nbsv);

Sigmaold  = Sigma;
while loop 
    nloop = nloop+1;
    %---------------------------------------------
    % Update Sigma
    %---------------------------------------------
    [Sigma,Alpsup,w0,pos,nbsv,obj0] = myMKLmulticlassupdate(K,Sigmaold,pos,Alpsup,w0,...
        C,yapp,nbsv,grad,obj0);
    %----------------------------------------------------
    % process approximate KKT conditions
    %----------------------------------------------------
    [grad] = myGradsvmoneagainstall(K,pos,Alpsup,nbsv);
    %------------------------------------------
    %  verbosity
    %------------------------------------------
    if verbose
        if nloop == 1 || rem(nloop,10)==0
            fprintf('--------------------------------------\n');
            fprintf('Iter | Obj.    | DiffBetas | KKT C    |\n');
            fprintf('---------------------------------------\n');
        end;
        fprintf('%d   | %8.4f | %6.4f   | %6.4f |\n',[nloop obj   max(abs(Sigma-Sigmaold)) KKTconstraint]);
    end
    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  max(abs(Sigma - Sigmaold))<1e-5 || nloop>=nloopmax
        loop = 0;
        fprintf(1,'variation convergence criteria reached \n');
    end
    %----------------------------------------------------
    % Updating Variables
    %----------------------------------------------------
    Sigmaold  = Sigma;
end