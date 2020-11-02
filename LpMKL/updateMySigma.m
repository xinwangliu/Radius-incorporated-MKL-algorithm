function [Sigma] = updateMySigma(w,nbsv,pos,K,Sigma0,qnorm)

numKer = size(K,3);
nbclass = length(nbsv);
nbsv=[0 nbsv];
aux=cumsum(nbsv);

f_tm2 = zeros(numKer,nbclass);
for ik = 1:numKer
    for ic = 1:nbclass
        %% Kernel matrix is given as a parameter
        indx   = aux(ic)+1:aux(ic)+nbsv(ic+1);
        waux   = w(indx);
        posaux = pos(indx);
        f_tm2(ik,ic) = Sigma0(ik)^2*waux'*K(posaux,posaux,ik)*waux;
    end
end
sum_f_tm2 = sum(f_tm2,2);
norm_ker_q = sum_f_tm2.^(1/(1+qnorm));
Sigma = norm_ker_q/sum(norm_ker_q.^qnorm).^(1/qnorm);