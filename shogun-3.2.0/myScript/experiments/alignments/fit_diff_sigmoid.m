%plots talignments and sigmoid fit
figure(7);
clf;
title('diff_unsorted_alignments');

min_y = 0;


y= talignments;
y=sort(y,'descend');	
y = abs(diff([y 0 ]));
x=[1:length(y)];
plot(y);
hold on;
title(sprintf('Scenario %d',i));
max_y = max(y);
%plot( data.no_inf(i)*ones(1,2), [min_y max_y], 'r', 'linewidth',2);
%f = @(p,x) p(1) + p(2) ./ (1 + exp((x-p(3))/p(4)));
%p = nlinfit(x,y,f,[0 0.02 5 1]);
%line(x,f(p,x),'color','g','linewidth',2);
a=axis;
axis([ a(1:2) min_y max_y]);
	