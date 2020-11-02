%plots talignments and sigmoid fit
figure(3);
clf;

min_y = 0;

yy = talignments(:)';
xx = [1:length(yy)];
plot(yy);
hold on;
max_y = max(yy);
%f = @(p,x) p(1) + p(2) ./ (1 + exp((xx-p(3))/p(4)));
%p = nlinfit(xx,yy,f,[0 0.2 50 10]);
%line(xx,f(p,xx),'color','g','linewidth',2);
%a=axis;
%axis([ a(1:2) min_y max_y]);
%plot( p(3)*ones(1,2), [min_y max_y], 'r', 'linewidth',2);
title('unsorted fit of talignments')

%fit(x,y,'a + b ./ (1 + exp(-(x-m)/s))','start',[0 20 50 500]) 