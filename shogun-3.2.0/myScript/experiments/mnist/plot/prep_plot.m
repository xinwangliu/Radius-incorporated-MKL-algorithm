function prep_plot(dat, col, idx, num,c)

set(gca,'Fontsize',21)
cmap=colormap;
c=cmap(ceil((idx-1)*size(cmap,1)/num+1), :);

col='-';
%if mod(idx,2)==0,
%	col='--';
%end

%cols={'k-d','m-d','c:o','c--x','c-.*','c-d','m-','b-d','b-v','b-^','b-<','->','-p','-h','-','--','-o','-x','-+','-*','-s','-d','-v','-^','-<','->','-p','-h'};
%col=cols{idx};
%if mod(idx,3)==2,
	%col=':';
%end

lw = 3;
%ms = 8.5;
[i,j,s]=find(dat);
x=unique(i);
y=zeros(1,length(x));
v=zeros(1,length(x));
for idx=1:length(x),
	d=nonzeros(dat(x(idx),:));
	y(idx)=mean(d);
	v(idx)=std(d);
	%keyboard
	idx=idx+1;
end
%keyboard
%y=mean(dat,2);
%v=std(dat');
%x=find(y);
%y=full(y(x));
%v=full(v(x));
%keyboard
%loglog(x,y,col,'LineWidth',lw,'MarkerSize',ms);
%loglog(x,y,col,'LineWidth',lw, 'Color', c,'MarkerSize',ms);
errorbar(x,y,v,col,'LineWidth',lw, 'Color',c);
errorbarlogx
errorbarlogy
hold on;
