% Merge many behavior mat files into a single one
clear all
cd('C:\Scratch\Beh');
%Get files
  [filename,pathname] = uigetfile(fullfile(pwd,'Data/Jean','*.mat'),'Select the sessions to append','MultiSelect','on');
c = length(filename);
%first iteration
    load(strcat(pathname,filename{1}))
    numtrials = length(trialData.trialEnd);
    tfields = fieldnames(trialData);
    for k=1:length(tfields) %clean trials that didnt happen
            if (isnumeric( trialData.(tfields{k}) ) || iscell( trialData.(tfields{k}) )) && length(trialData.(tfields{k}))>numtrials
                trialData.(tfields{k}) = trialData.(tfields{k})(1:numtrials);
            end
    end
    trialDataM = trialData;
    MergedStr{1} = filename{1};
    clear trialData sessionData
    
%
if c>1
    for k = 2:c
      load(strcat(pathname,filename{k}))
      MergedStr{k} = filename{k};
      numtrials = length(trialData.trialEnd);
      tfields = fieldnames(trialData);
      for kk=1:length(tfields) %Build new trialDataM
            if isnumeric( trialData.(tfields{kk}) ) || iscell( trialData.(tfields{kk}) )
                if length(trialData.(tfields{kk}))>numtrials %^clean trials that didnt happen
                    trialData.(tfields{kk}) = trialData.(tfields{kk})(1:numtrials);
                end    
                trialDataM.(tfields{kk}) = [ trialDataM.(tfields{kk}) ; trialData.(tfields{kk})];
            end
      end
    end
end
trialData = trialDataM;
clear trialDataM
%Save the new variables into a new merged file
namenew = inputdlg('Enter name of created merged file');
cd('C:\Scratch\Beh\Data')
save(namenew{1},'sessionData','trialData','MergedStr')
clear all