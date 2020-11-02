figs = findobj('type','figure');
for j=1:length(figs)
  i= figs(j);
  saveas(i, sprintf('~/temp/plot%d',i),'fig');
  saveas(i, sprintf('~/temp/plot%d',i),'eps');
  saveas(i, sprintf('~/temp/plot%d',i),'pdf');
end
