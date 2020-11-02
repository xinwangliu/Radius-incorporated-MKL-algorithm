function KC = centerlizedKernel(K)
[n1,n2]= size(K);

e1 = ones(n1,1);
em1 = e1*e1'/n1;
in1 = eye(n1);
cm1 = (in1 - em1);

e2 = ones(n2,1);
em2 = e2*e2'/n2;
in2 = eye(n2);
cm2 = (in2 - em2);
KC = cm1*K*cm2;