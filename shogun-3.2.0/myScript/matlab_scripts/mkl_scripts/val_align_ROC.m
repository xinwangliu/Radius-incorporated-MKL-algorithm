function tps = val_align_ROC(roc, fps)
% val_align_ROC - compute an "aligned" ROC-curve over several trials
%
% Synopsis:
%   tps = val_align_ROC(roc, fps)
%
% Arguments:
%   roc:  a ROC-curve structure array ...
%         .tp   true positive vector
%         .fp   false positive vector
%   fps:  a vector of false-positive "sample points"
%
% Returns:
%   tps:  a matrix of true positive values at "sample points"
%
% Description: 
%   val_align_ROC aligns several ROC curves, possibly of different length,
%   obtained from different experiments. The ROC curves are assumed to be
%   stored as structure array with fields .tp and .fp. Alignment is
%   performed for a vector of "sample points" providing the false positive
%   values, for which the true positive values are searched for in all ROC
%   curves. The result is a [length(roc), length(fps)] matrix containing in
%   colums true positives from all experiment for a given false-positive
%   value.
%
% See also:
%   val_ROC
%
% $Id: val_align_ROC.m 13220 2007-10-16 14:56:40Z neuro_cvs $
%
% Copyright (C) 2005,2007 Fraunhofer FIRST
% Author: Konrad Rieck (rieck@first.fhg.de)
%         Pavel Laskov (laskov@first.fhg.de)

tps = [];

for a=1:length(roc)
  tpa = roc(a).tp;
  fpa = roc(a).fp;
  for t=1:length(fps)

    j = find(fpa >= fps(t));
    if ~isempty(j)
      [b, max_idx] = min(fpa(j));
      max_idx = j(max_idx);
    else
      [m,max_idx] = max(fpa);
    end

    j = find(fpa <= fps(t));
    if ~isempty(j)
       [b, min_idx] = max(fpa(j));
       min_idx = j(min_idx);
    else
      [m,min_idx] = min(fpa);
    end

    if fpa(max_idx) ~= fpa(min_idx)
       dy = tpa(max_idx) - tpa(min_idx);
       dx = fpa(max_idx) - fpa(min_idx);
       s = dy / dx;
       tps(a,t) = s * (fps(t) - fpa(min_idx)) + tpa(min_idx);
    else
       tps(a,t) = tpa(max_idx);
    end  

% Over emphasizes ROC area!
%    next_idx = max(find(fpa >= fps(t)));
%    [m,ind_max] = max(fpa);
%    if isempty(next_idx)
%      tps(a,t) = tpa(ind_max);
%    else
%      tps(a,t) = tpa(next_idx);
%    end
  end	
end;   
