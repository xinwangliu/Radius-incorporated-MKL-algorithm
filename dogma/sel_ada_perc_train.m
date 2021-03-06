function model = sel_ada_perc_train(X,Y,model)
% SEL_ADA_PERC_TRAIN Kernel Selective Perceptron algorithm, with adaptive
% parameter
%
%    MODEL = SEL_ADA_PERC_TRAIN(X,Y,MODEL) trains an classifier according
%    to the Selective Perceptron algorithm.
%
%    Additional parameters: 
%    - model.bs governs the sampling rate of the algorithm.
%      Default value is 1.
%
%   References:
%     - Cesa-Bianchi, N., Gentile, C., & Zaniboni, L. (2006).
%       Worst-Case Analysis of Selective Sampling for Linear Classification
%       Journal of Machine Learning Research, 7, (pp. 1205-1230).

%    This file is part of the DOGMA library for MATLAB.
%    Copyright (C) 2009-2011, Francesco Orabona
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    Contact the author: francesco [at] orabona.com

n = length(Y);   % number of training samples

if isfield(model,'iter')==0
    model.iter=0;
    model.w=zeros(1,size(X,1));
    model.w2=zeros(1,size(X,1));
    model.errTot=0;
    model.numSV=zeros(numel(Y),1);
    model.aer=zeros(numel(Y),1);
    model.pred=zeros(numel(Y),1);
    
    model.numQueries=0;
    model.maxR2=0;
    model.numUpdates=0;
end

if isfield(model,'bs')==0
    model.bs=1;
end

for i=1:n
    model.iter=model.iter+1;
        
    val_f=model.w*X(:,i);

    Yi=Y(i);
    
    model.errTot=model.errTot+(sign(val_f)~=Yi);
    model.aer(model.iter)=model.errTot/model.iter;
    model.pred(model.iter)=val_f;
    
    R2=max(norm(X(:,i))^2,model.maxR2);
    b=model.bs*R2*sqrt(1+model.numUpdates);
    Z=(rand<b/(abs(val_f)+b));
    
    model.numQueries=model.numQueries+Z;
    
    if Z==1 && Yi*val_f<=0
        model.w=model.w+Yi*X(:,i)';
        model.S(end+1)=model.iter;
        model.numUpdates=model.numUpdates+1;
        model.maxR2=R2;
    end

    model.w2=model.w2+model.w;
    
    model.numSV(model.iter)=numel(model.S);
    
    if mod(i,model.step)==0
      fprintf('#%.0f SV:%5.2f(%d)\tQueried Labels:%5.2f(%d)\tAER:%5.2f\n', ...
          ceil(i/1000),numel(model.S)/model.iter*100,numel(model.S),...
          model.numQueries/model.iter*100,model.numQueries,...
          model.aer(model.iter)*100);
    end
end
