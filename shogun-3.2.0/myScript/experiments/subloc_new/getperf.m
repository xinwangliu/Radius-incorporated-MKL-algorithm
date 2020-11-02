function [confmat,accuracy,tp,fp,fn,tn,precision,recall,F,MCC] = getperf(actual,pred,classes)
% getperf : gets confusion matrices, precision, recall, F scores,
%         true positives, false positives, false negatives
% [confmat,accuracy,tp,fp,fn,tn,precision,recall,F,MCC] = getperf (actual,pred,[classes])
%
% actual is a N-element vector representing the actual classes
% pred is a N-element vector representing the predicted classes
% classes is a vector with the numbers of the classes (by default, it is 0:k, where k is the
%    largest integer to appear in actual or pred.
%

% All programs in this collection are free software: 
% you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright 2007 Alexander Zien and Cheng Soon Ong


if size(actual,1) ~= size(pred,1)
  pred=pred';
end
if nargin < 3
  classes = [0:max(max(actual),max(pred))];
end

m = length(actual);

%accuracy = sum(actual==pred)/m;
for i=1:length(classes)
  a = classes(i);
  d = find(actual==a);     % d has indices of points with class a
  for j=1:length(classes)
    confmat(i,j) = length(find(pred(d)==classes(j)));
  end
end

accuracy=[];
precision=[];
recall=[];
F=[];
tp = [];
fp = [];
fn = [];
tn = [];
MCC = [];
for i=1:length(classes)
  S =  sum(confmat(i,:));
  if S > 0
    recall(i) = confmat(i,i) / sum(confmat(i,:));
  else
    recall(i) = 0;
  end
  tp(i) = confmat(i,i);
  S =  sum(confmat(:,i));
  fp(i) = S - confmat(i,i);
  fn(i) = sum(confmat(i,:),2) - confmat(i,i);
  d = diag(confmat);
  d(i) = 0;
  tn(i) = sum(d);
  if S
    precision(i) = confmat(i,i) / sum(confmat(:,i));
  else
    precision(i) = 0;
  end
  if (precision(i)+recall(i))
    F(i) = 2 * (precision(i)*recall(i)) / (precision(i)+recall(i));
  else
    F(i) = 0;
  end
  
  S = tp(i)+tn(i)+fp(i)+fn(i);
  if S
    accuracy(i) = (tp(i)+tn(i))/S;
  else
    accuracy(i) = 0;
  end
  
  
  S = (tp(i)+fn(i))*(tp(i)+fp(i))*(tn(i)+fp(i))*(tn(i)+fn(i));
  if S
    MCC(i) = (tp(i)*tn(i)-fp(i)*fn(i))/sqrt(S);
  else
    MCC(i) = 0;
  end
  
end


%resmat = [tp',fp',fn',precision',recall',F'];
