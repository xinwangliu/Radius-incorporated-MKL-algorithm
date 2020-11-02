
clear;

datasetNames = {};
datasetNames{end+1} = 'plant';
%datasetNames{end+1} = 'nonpl';
%datasetNames{end+1} = 'psortNeg';
%datasetNames{end+1} = 'psortPos';

%cplexLicence = '587380';
cplexLicence = '820980';

methodNames = {};
methodNames{end+1} = [ 'CPLEX-' cplexLicence ];
methodNames{end+1} = 'NEWTON1';
methodNames{end+1} = 'DIRECT';

%epsSettings = { 'b' };
epsSettings = { 'a', 'b' };

nofSplits = 30;



for( i = 1:length(datasetNames) )
  datasetName = datasetNames{ i };
  for( j = 1:length(methodNames) )
    methodName = methodNames{ j };
    for( k = 1:length(epsSettings) )
      epsSetting = epsSettings{ k };
      % - eval
      evalResult;
      %close all;
    end;
  end;
end;



