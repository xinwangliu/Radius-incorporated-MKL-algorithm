function [Sigma,Alpsup,w0,pos,nbsv] = myMKLmulticlass(K,yapp,C,verbose)


nbkernel=size(K,3);
nloopmax=1000;
option.algo='oneagainstall';
option.verbosesvm=0;
Sigma=ones(nbkernel,1)/nbkernel;
%--------------------------------------------------------------------------------
% Options used in subroutines
%--------------------------------------------------------------------------------
option.goldensearch_deltmax=1e-1;
option.goldensearchmax=1e-8;
option.firstbasevariable='first';
option.numericalprecision=1e-8;

%------------------------------------------------------------------------------%
% Initialize
%------------------------------------------------------------------------------%
nloop = 0;
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
        C,yapp,nbsv,grad,obj0,option);
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
            fprintf('Iter | Obj.    | DiffBetas |\n');
            fprintf('---------------------------------------\n');
        end;
        fprintf('%d   | %8.4f | %6.4f |\n',[nloop obj0   max(abs(Sigma-Sigmaold))]);
    end
    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  max(abs(Sigma - Sigmaold))<1e-5 || nloop>=nloopmax
        loop = 0;
        fprintf(1,'variation convergence criteria reached\n');
    end
    %----------------------------------------------------
    % Updating Variables
    %----------------------------------------------------
    Sigmaold  = Sigma;
end