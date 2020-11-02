function [Weigth,InfoKernel]=UnitTraceNormalization(x,kernelvec,kerneloptionvec,variablevec)

%%%
display('Kernel Normalization')
N=size(x,1);  %the number of samples
chunksize=200;
nbk=1;
for i=1:length(kernelvec);
    % i
    for k=1:length(kerneloptionvec{i})

        somme=0;

        chunks1=ceil(N/chunksize);

        for ch1=1:chunks1
            ind1=(1+(ch1-1)*chunksize) : min( N, ch1*chunksize);
            if(0)
             tmpK = svmkernel(x(ind1,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k));
             somme = somme + (sum(diag(tmpK))- sum(sum(tmpK))/size(tmpK,1));
            else
              somme=somme+sum(diag(svmkernel(x(ind1,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k))));
            end
        end;
        %(trace(tmpK)) %- sum(sum(tmpK))/size(tmpK,1)
        %         for j=1:N
        %             somme=somme+svmkernel(x(j,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k));
        %
        %         end
        if somme~=0
            Weigth(nbk)=1/somme;
            InfoKernel(nbk).kernel=kernelvec{i};
            InfoKernel(nbk).kerneloption=kerneloptionvec{i}(k);
            InfoKernel(nbk).variable=variablevec{i};
            InfoKernel(nbk).Weigth=1/somme;
            nbk=nbk+1;
%         else
%             A
        end;
    end;
end;