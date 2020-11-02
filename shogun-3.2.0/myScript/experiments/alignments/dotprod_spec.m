function p=dotprod_spec(K,M)

p= max(eig(K*M));

%p = trace(K*M);  %frob