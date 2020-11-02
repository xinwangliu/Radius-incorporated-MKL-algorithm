function [tp,fp,auc] = val_ROC(s,Ytr,varargin)
% val_ROC - compute a receiver operating characteristic (ROC)
%
% Synopsis: 
%   [tp,fp,auc] = val_ROC(s,Ytr)
%   [tp,fp,auc] = val_ROC(s,Ytr,'property',value)
% 
% Arguments: 
%   s:   anomaly scores 
%   Ytr: labels (>=1 for attacks, 0 for normal)
%
% Returns:
%    tp:  an array of true-positive values
%    fp:  an array of false-positive values
%    auc: area under curve
% 
% Properties:
%    'ub':   Specify upper bound for the false-positive rate
%    'lb':   Specify lower bound for the false-positive rate
%  
% $Id: val_ROC.m 9786 2006-01-02 09:30:22Z neuro_cvs $
% 
% Copyright (C) 2004 Fraunhofer FIRST
% Author: Pavel Laskov (laskov@first.fhg.de), 
%         Konrad Rieck (rieck@first.fhg.de)
%         Patrick Duessel (duessel@first.fhg.de)

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, 'lb',0, 'ub', 1.0);

% Simplify labels
Ytr(find(Ytr >= 1)) = 1;
Ytr(find(Ytr <= 0)) = 0;

% Determine indices
indp = find(Ytr==1);
indm = find(Ytr==0);

T = sort(s(indp));
k = length(T);

tp = zeros(1,k+2);
fp = zeros(1,k+2);

tp(1) = 1;
fp(1) = 1;
for j = 1:k
  indc = find(s>=T(j));
  tp(j+1) = sum(Ytr(indc))/length(indp);
  fp(j+1) = sum(1-Ytr(indc))/length(indm);
end
tp(end) = 0;
fp(end) = 0;

fidx_ub = find(fp < opt.ub);
fidx_lb = find(fp > opt.lb);
fidx = intersect(fidx_ub, fidx_lb);

%%% find previous and next element beyond bounds %%%
i1 = max(find(fp >= opt.ub));
i2 = min(find(fp <= opt.lb));

ub_nfp = fp(i1); 
ub_ntp = tp(i1); 
lb_nfp = fp(i2);
lb_ntp = tp(i2);

range = opt.ub - opt.lb;

tp = tp(fidx);
fp = fp(fidx);

%%% add first and last approximated value to fp, tp

if ~isempty(fidx)          
  lb.tp = [tp(end) lb_ntp];
  lb.fp = [fp(end) lb_nfp];

  ub.tp = [tp(1) ub_ntp];
  ub.fp = [fp(1) ub_nfp];
else
  % there is no data point within interval bound
  lb.tp = [ub_ntp lb_ntp];
  lb.fp = [ub_nfp lb_nfp];

  ub.tp = [lb_ntp ub_ntp];
  ub.fp = [lb_nfp ub_nfp];
end

[ftp, ltp] = approx_ROC(lb, ub, opt);
tp(2:end+1) = tp;	
fp(2:end+1) = fp;

tp(1) = ltp;
fp(1) = opt.ub;
tp(end+1) = ftp;
fp(end+1) = opt.lb;

auc = -diff(fp)*tp(2:end)' / range;

%%% returns first and last tp additionally added to tp list
function [ftp, ltp] = approx_ROC(lb, ub, opt)

  %%% upper bound approximation %%%
  m = (ub.tp(2) - ub.tp(1)) / (ub.fp(2) - ub.fp(1));
  ltp = m*opt.ub + ub.tp(1) - m*ub.fp(1);

  %%% lower bound approximation %%%
  m = (lb.tp(1) - lb.tp(2)) / (lb.fp(1) - lb.fp(2));
  ftp = m*opt.lb + lb.tp(2) - m*lb.fp(2); 
