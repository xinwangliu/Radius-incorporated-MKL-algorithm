figure(5);
clf;

for  i=1:5
  subplot(2,3,i);
  e = eig(K(:,:,i));
  e=flipud(e(:));
  eigens(:,i) = e;
  plot(e);
end  