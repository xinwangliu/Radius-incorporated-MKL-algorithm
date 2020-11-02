function [prc_auc, roc_auc, tp, fp, ppv] = val_aucs(s,Ytr,varargin)
% val_ROC - compute a positive preditive value and ROC AUC
%
% Synopsis: 
%   [prc_auc,roc_auc,tp,fp,ppv] = val_aucs(s,Ytr)
%   [prc_auc,roc_auc,tp,fp,ppv] = val_aucs(s,Ytr,'property',value)
% 
% Arguments: 
%   s:   anomaly scores 
%   Ytr: labels (>=1 for attacks, 0 for normal)
%
% Returns:
%    ppv: an array of ppv values
%    fp:  an array of false-positive values
%    auc: area under curve
% 
% $Id: val_aucs.m 14630 2008-10-21 17:18:02Z neuro_cvs $
% 
% Copyright (C) 2006 Fraunhofer FIRST
% Author: Alex Ziehn
%         Konrad Rieck (rieck@first.fhg.de)

opt = propertylist2struct(varargin{:});
opt = set_defaults(opt, 'gra', 10000);

if size(Ytr,1) > size(Ytr,2)
   Ytr = Ytr';
end

if size(s,1) > size(s,2)
   s = s';
end

% Simplify labels
Ytr(find(Ytr >= 1)) = 1;
Ytr(find(Ytr <= 0)) = 0;

% Count label values
n = length(Ytr);
np = sum(Ytr == 1);
nn = sum(Ytr == 0);

% Outpos 
outPos = s( Ytr == 1 );
outNeg = s( Ytr == 0 );


% Sort scores
[a perm] = sort(s);

TP = np;
TN = 0;
FP = nn;
FN = 0;
TPRs = nan * zeros( n+1, 1 );
FPRs = nan * zeros( n+1, 1 );
PPVs = nan * zeros( n+1, 1 );
TPRs(1) = 1;
FPRs(1) = 1;
PPVs(1) = 0;

for i = 1:n
   j = perm(i);
   if Ytr(j) == 0
     FP = FP - 1;
     TN = TN + 1;
   else
     TP = TP - 1;
     FN = FN + 1;
   end  
   sn = TP / ( TP + FN );
   sp = TN / ( TN + FP );
   TPRs(i+1) = sn;
   FPRs(i+1) = 1 - sp;
   if( TP+FP == 0 )
      PPVs(i+1) = 1;
   else
      PPVs(i+1) = TP / ( TP + FP );
   end
end

granularity = opt.gra;

% === compute ROC
thresholds = -sort( -outPos );
if( np < 2*granularity )
  p = np;
else
  p = granularity;
end;
TPRs = (1:p) / p;
FPRs = nan*zeros( p, 1 );
PPVs = nan*zeros( p, 1 );
for( i = 1:p )
  tpr = TPRs( i );
  r = ceil(np/p*i);
  if r > length(thresholds)
     r = length(thresholds);
  end
  if r < 1 
     r = 1;
  end
  threshold = thresholds(r);
  tp = sum( outPos >= threshold );
  fp = sum( outNeg >= threshold );

  % --- correct for ties
  nofPosMatch = sum( outPos == threshold );
  nofNegMatch = sum( outNeg == threshold );
  TP = np * tpr;
  gapPos = round( tp - TP );
  if( gapPos > 0 )
    frac = gapPos / nofPosMatch;
    gapNeg = frac * nofNegMatch;
    tp = TP;
    fp = fp - gapNeg;
  end;
  % --- store point
  FPRs(i) = fp / nn;
  PPVs(i) = tp / ( tp + fp );
end;

roc_auc = 1 - mean(FPRs);
prc_auc = mean(PPVs);
tp = TPRs;
fp = FPRs;
ppv = PPVs;
