function targetmat = construct_targetmat(labels1,labels2)

   n = length(labels1);
   m = length(labels2);
   targetmat_same = zeros(n,n);
   targetmat_diff = zeros(n,n);
   for i=1:n
      for j=1:m
	 if (labels1(i)==labels2(j))
	    targetmat_same(i,j) = 1;
	 end
	 if (labels1(i)~=labels2(j))
	    targetmat_diff(i,j) = 1;
	 end
      end
   end
   targetmat = targetmat_same - targetmat_diff; % + eye(n);

 
