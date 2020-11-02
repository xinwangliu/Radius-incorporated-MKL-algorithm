function [Weigth,InfoKernel]=UnitTraceNormalization2(x,kernelvec,kerneloptionvec,variablevec)

%aa = [0.7136    0.5531    0.5524    0.6250    0.7008    0.1920    0.9114    0.7452    0.6574    0.5269];
aa = [0.3336    0.1109    0.7096    0.8590    0.6025    0.4741    0.1283    0.8514    0.9937    0.0888];
%aa = [0.5381    0.4757    0.6979    0.8259    0.1634    0.7648    0.6213    0.6983    0.5411    0.9443];
%aa =[0.1173    0.8353    0.9276    0.3768    0.8073    0.5174    0.0571    0.3642    0.1486    0.3888];
%aa =[0.0733    0.1454    0.5832    0.8124    0.2774    0.0428    0.1359    0.8705    0.8602    0.2956];

% aa = [0.7330    0.7949    0.1445    0.8572    0.6943    0.4358    0.0781    0.4036    0.3973    0.2942];
% aa = [0.0801    0.0546    0.7302    0.8572    0.8968    0.6004    0.6245    0.4803    0.5700    0.2635];
% aa = [0.1853    0.7038    0.9609    0.8268    0.2441    0.8896    0.4626    0.0236    0.4816    0.9785];
% aa = [0.5065    0.8948    0.8964    0.2869    0.0403    0.4980    0.2220    0.2426    0.8381    0.3760];
% aa = [0.8021    0.6773    0.5369    0.6065    0.9588    0.3655    0.7388    0.8150    0.1655    0.7032];
% aa = [0.8283    0.45031    0.9125    0.3926    0.8705    0.9743    0.1834    0.3602    0.6658    0.2810];

N=size(x,1);
chunksize=N;
nbk=1;
for i=1:length(kernelvec);
    % i
    for k=1:length(kerneloptionvec{i})

        somme=0;

        chunks1=ceil(N/chunksize);

        for ch1=1:chunks1
            ind1=(1+(ch1-1)*chunksize) : min( N, ch1*chunksize);
            somme=somme+sum(diag(svmkernel(x(ind1,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k))));
        end;
        %         for j=1:N
        %             somme=somme+svmkernel(x(j,variablevec{i}),kernelvec{i},kerneloptionvec{i}(k));
        %
        %         end
        if somme~=0
            a=aa(nbk); 
            Weigth(nbk)=a;
            InfoKernel(nbk).kernel=kernelvec{i};
            InfoKernel(nbk).kerneloption=kerneloptionvec{i}(k);
            InfoKernel(nbk).variable=variablevec{i};
            InfoKernel(nbk).Weigth=a;
            nbk=nbk+1
%         else
%             A
        end;
    end;
end;