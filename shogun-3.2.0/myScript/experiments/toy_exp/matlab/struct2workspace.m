function struct2workspace(s);
% takes deep contents of structure, and
% floats to top of workspace
% Scott Seidman, 7/05
names=fieldnames(s);
for i=1:length(names)
assignin('caller',names{i},s.(names{i}));
end
