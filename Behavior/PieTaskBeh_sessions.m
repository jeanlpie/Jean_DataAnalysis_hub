%% PieTaskBeh_sessions
%Select n sessions. This code runs behavioral analyses with statistics
%across sessions.

%% Import files
clear all
cd('C:\Scratch\Beh');
%Get files
  [filename,pathname] = uigetfile(fullfile(pwd,'Data/Jean','*.mat'),'Select the sessions to append','MultiSelect','on');
numses = length(filename);

%% First overview off all sessions, to cut out trials before, in-the-middle or at the end
DataSessions = struct;

for ises = 1:numses
    load(strcat(pathname,filename{ises}))
    numtrials = length(trialData.correctResponse);
    %Show performance
    sliding_perf_fun(sessionData,trialData,40);
    % Show start of session
    %view_it(trialData,sessionData,1,20,[-1.500 1.5]);
    cut_start(ises) = input(strcat('Enter number of trials to trim at START of session',filename{ises}));
    % Show end of session
    %view_it(trialData,sessionData,numtrials-30,numtrials,[-1.500 1.5]);
    cut_end(ises) = input(strcat('Enter trial number to cut session at the END','- numtrials:',num2str(numtrials),filename{ises}));
    % Cut trials
    tfields = fieldnames(trialData);
    for k=1:length(tfields) %clean trials that didnt happen
            if (isnumeric( trialData.(tfields{k}) ) || iscell( trialData.(tfields{k}) )) && length(trialData.(tfields{k}))>numtrials
                trialData.(tfields{k}) = trialData.(tfields{k})(1+cut_start(ises):cut_end(ises));
            end
    end
    close all
    DataSessions(ises).trialData = trialData;
    DataSessions(ises).sessionData = sessionData;
end
%Save the new structure
namenew = inputdlg('Enter name of created merged file');
cd('C:\Scratch\Beh\Data')
save(namenew{1},'DataSessions')
clear all

%% Now the real analysis
%Prelim
onlyvis= [7 6 5 1 2 3 4];       onlyvis_opto = [16 17 18];
onlytac = [11 10 1 8 9];        onlytac_opto = [16 19];
facilvis = [15 14 1 12 13];     facilvis_opto = [16 20 21];

numses = size(DataSessions,2);
R = struct;

SelSes = [3 4 6 7 8 9];%[1 3:9]; %1:numses %[1 3 4 5 6 7]; %1:numses
%Loop through sessions

for ises = SelSes
trialData = DataSessions(ises).trialData;
sessionData = DataSessions(ises).sessionData;    
%Determine if it's first or second version
if isfield(trialData,'conditionnumber')
    booliscond = true;
else 
    booliscond = false;
end
%Convert nans into zeros:
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
%Contrasts in log ?
% trialData.stimContrast = log10(trialData.stimContrast*1000);
% trialData.stimContrast(isinf(trialData.stimContrast)) = 0; %ATTENTION ! Contrast cannot be smaller or equal to 1.
%Make left negative
if isfield(trialData,'stimVISSide')
    trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L') =  trialData.stimDeflection(cell2mat(trialData.stimTACSide)=='L')*-1;
    trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')   =  trialData.stimContrast(cell2mat(trialData.stimVISSide)=='L')*-1;    
else
    trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
    trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
end
%Compute response ratios per cond

if booliscond
    
    conds = unique(trialData.conditionnumber);
    condsnummax = max(conds);
    for k = 1:condsnummax
          R(ises).TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
          R(ises).HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / R(ises).TotalTrials(k);
          R(ises).HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / R(ises).TotalTrials(k);
          R(ises).NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / R(ises).TotalTrials(k);
          R(ises).Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
          R(ises).Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
                %for fisher
          R(ises).RR(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k));
          R(ises).LL(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k));
          R(ises).NN(k) = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k));
    end 
else
    R(ises).Deflections = unique(trialData.stimDeflection); R(ises).Deflections(isnan(R(ises).Deflections))=[];
    R(ises).Contrasts = unique(trialData.stimContrast); R(ises).Contrasts(isnan(R(ises).Contrasts))=[];
    %  Compute percent right, left, and no response
    NTac = length(R(ises).Deflections);
    NVis = length(R(ises).Contrasts);
    for v=1:NVis
        for t=1:NTac
          R(ises).TotalTrials(v,t) =  nansum(trialData.stimDeflection==R(ises).Deflections(t) & trialData.stimContrast==R(ises).Contrasts(v));  
          R(ises).HitRight(v,t) = nansum(trialData.firstRespLFR=='R' & trialData.stimDeflection==R(ises).Deflections(t) & trialData.stimContrast==R(ises).Contrasts(v)) / R(ises).TotalTrials(v,t);
          R(ises).HitLeft(v,t)  = nansum(trialData.firstRespLFR=='L' & trialData.stimDeflection==R(ises).Deflections(t) & trialData.stimContrast==R(ises).Contrasts(v)) / R(ises).TotalTrials(v,t);
          R(ises).NoResp(v,t)   = nansum(isnan(trialData.firstRespLFR) & trialData.stimDeflection==R(ises).Deflections(t) & trialData.stimContrast==R(ises).Contrasts(v)) / R(ises).TotalTrials(v,t);
        end
    end
end

end
%% Fisher test
%fisher test
comp = [1,16 ; 2,17 ; 3,18 ; 8,19 ; 12,20 ; 13,21];
for k = 1:size(comp,1)
xL = [LL(comp(k,1)),TotalTrials(comp(k,1))-LL(comp(k,1)) ; LL(comp(k,2)) , TotalTrials(comp(k,2))- LL(comp(k,2)) ];
xR = [RR(comp(k,1)),TotalTrials(comp(k,1))-RR(comp(k,1)) ; RR(comp(k,2)) , TotalTrials(comp(k,2))- RR(comp(k,2)) ];
xN = [NN(comp(k,1)),TotalTrials(comp(k,1))-NN(comp(k,1)) ; NN(comp(k,2)) , TotalTrials(comp(k,2))- NN(comp(k,2)) ];
[hL(k),pL(k)] = fishertest(xL);
[hR(k),pR(k)] = fishertest(xR);
[hN(k),pN(k)] = fishertest(xN);
end

%% Plot Only VIS
figure()
hold on
for ises = SelSes %1:numses
    if isfield(trialData,'conditionnumber')
plot(R(ises).Contrasts(onlyvis), R(ises).HitRight(onlyvis),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(R(ises).Contrasts(onlyvis), R(ises).HitLeft(onlyvis),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(R(ises).Contrasts(onlyvis), R(ises).NoResp(onlyvis),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
    else
plot(R(ises).Contrasts,R(ises).HitRight(:,(R(ises).Deflections==0)),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(R(ises).Contrasts,R(ises).HitLeft(:,(R(ises).Deflections==0)),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(R(ises).Contrasts,R(ises).NoResp(:,(R(ises).Deflections==0)),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
    end
end
legend({'LickRight','LickLeft','NoResponse'})
ylim([0 1])
title('Mouse choices - Visual Only','FontSize',18)
xlabel('<-LEFT      log(10*Contrast)      RIGHT->','FontSize',18)
ylabel('Percent choice','FontSize',18)
set(gca,'FontSize',14)

%% Plot Only TAC
figure()
hold on
for ises = SelSes
        if isfield(trialData,'conditionnumber')
plot(R(ises).Deflections(onlytac), R(ises).HitRight(onlytac),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(R(ises).Deflections(onlytac), R(ises).HitLeft(onlytac),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(R(ises).Deflections(onlytac), R(ises).NoResp(onlytac),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
        else
plot(R(ises).Deflections,R(ises).HitRight((R(ises).Contrasts==0),:),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(R(ises).Deflections,R(ises).HitLeft((R(ises).Contrasts==0),:),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(R(ises).Deflections,R(ises).NoResp((R(ises).Contrasts==0),:),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)            
        end
        
end
legend({'LickRight','LickLeft','NoResponse'})
ylim([0 1])
title('Mouse choices - Tactile Only','FontSize',18)
ylabel('Percent choice','FontSize',18)
xlabel('<-LEFT      Deflection      RIGHT->','FontSize',18)
legend({'LickRight','LickLeft','NoResponse'})
set(gca,'FontSize',14)

%% IMAGESCS
figure('Units','normalized','pos',[0.05 0.4 0.8 0.5])
colormap(summer)
% for ises = SelSes
%     
% 
% end

subplot(1,3,1)
imagesc(HitLeft*100,'AlphaData',~isnan(HitLeft))
set(gca,'clim',[0,100],'XTick',1:length(Deflections),'XTickLabels',Deflections,'YTick',1:length(Contrasts),'YTickLabels',Contrasts);
xlabel('TAC')
ylabel('VIS')
title('Percent choice Left')
posx=[];
for tt=1:NTac
posx = [ posx repmat(tt,1,NVis)];
end
posy = [repmat(1:NVis,1,NTac)];
text(posx(~isnan(HitLeft)),posy(~isnan(HitLeft)),num2str((round(( HitLeft(~isnan(HitLeft(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')



subplot(1,3,2)
imagesc(HitRight*100,'AlphaData',~isnan(HitRight))
set(gca,'clim',[0,100],'XTick',1:length(Deflections),'XTickLabels',Deflections,'YTick',1:length(Contrasts),'YTickLabels',Contrasts);
xlabel('TAC')
ylabel('VIS')
title('Percent choice Right')
text(posx(~isnan(HitRight)),posy(~isnan(HitRight)),num2str((round(( HitRight(~isnan(HitRight(:))) ),2)*100)),'Color',[0.3 0.1 0],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')



subplot(1,3,3)
imagesc(NoResp*100,'AlphaData',~isnan(NoResp))
set(gca,'clim',[0,100],'XTick',1:length(Deflections),'XTickLabels',Deflections,'YTick',1:length(Contrasts),'YTickLabels',Contrasts);
xlabel('TAC')
ylabel('VIS')
title('Percent NoResponse')
text(posx(~isnan(NoResp)),posy(~isnan(NoResp)),num2str((round(( NoResp(~isnan(NoResp(:))) ),2)*100)),'Color',[.15 .15 .15],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')

cb = colorbar('Units','normalized','Position',[0.92 0.1 0.02 0.8]);
cb.Label.String = 'Choice percentage';


%% Make it one long session for other stuff
%load DataSessions again

c = size(DataSessions,2);
%first iteration
trialData = DataSessions(1).trialData;
sessionData = DataSessions(1).sessionData;
    numtrials = length(trialData.trialEnd);
    tfields = fieldnames(trialData);
    trialDataM = trialData;
    clear trialData sessionData
    
SelSes = 1:c; %[3 4 6 7 8 9] %[1 3:9] %[1 3:7]; %1:c    
    
    for ises = SelSes
        trialData = DataSessions(ises).trialData;
        sessionData = DataSessions(ises).sessionData;
      numtrials = length(trialData.trialEnd);
      tfields = fieldnames(trialData);
      for k=1:length(tfields) %Build new trialDataM
            if isnumeric( trialData.(tfields{k}) ) || iscell( trialData.(tfields{k}) )
                trialDataM.(tfields{k}) = [ trialDataM.(tfields{k}) ; trialData.(tfields{k})];
            end
      end
    end

trialData = trialDataM;
clear trialDataM c k 
%% Preproc
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;
% %Contrasts in log ?
trialData.stimContrast = log10(trialData.stimContrast*1000);
trialData.stimContrast(isinf(trialData.stimContrast)) = 0; %ATTENTION ! Contrast cannot be smaller or equal to 1.
%Make left negative
trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;
%Supposing M,V and T trials had the same intensity values:
Deflections = unique(trialData.stimDeflection); Deflections(isnan(Deflections))=[];
Contrasts = unique(trialData.stimContrast); Contrasts(isnan(Contrasts))=[];
%  Compute percent right, left, and no response
NTac = length(Deflections);
NVis = length(Contrasts);
for v=1:NVis
    for t=1:NTac
      TotalTrials(v,t) =  nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));  
      HitRight(v,t) = nansum(trialData.firstRespLFR=='R' & trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v)) / nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));
      HitLeft(v,t)  = nansum(trialData.firstRespLFR=='L' & trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v)) / nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));
      NoResp(v,t)   = nansum(isnan(trialData.firstRespLFR) & trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v)) / nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));
    end
end
%% Fit sigmoids
ft2 = fittype('g+(1-g-l)*(1/2)*(1+erf((x-m)/(s*sqrt(2))))');
coeffnames(ft2); % g l m s
%VIS
x = Contrasts;
yR = HitRight(:,(Deflections==0));
yL = HitLeft(:,(Deflections==0));
startPointsR = [ 0.15 , 0.2 , 1.9, 1];
startPointsL = [ 0.15 , 0.2 , -1.9, 1];
VWeights = TotalTrials(:,(Deflections==0));
fvisR = fit(x,yR,ft2,'Start', startPointsR, 'Weights', VWeights);
fvisL = fit(x,yL,ft2,'Start', startPointsL, 'Weights', VWeights);
VisCoeffsR = coeffvalues(fvisR);
VisCoeffsL = coeffvalues(fvisL);
%plot
figure();
hRv = plot(fvisR,x,yR);
set(hRv,'Color',[.8 0 0],'LineWidth',2)
hold on
hLv = plot(fvisL,x,yL);
set(hLv,'Color',[0 .8 0],'LineWidth',2)
ylim([0 1])
scatter(Contrasts,HitRight(:,(Deflections==0)),VWeights,[.8 0 0],'filled')
scatter(Contrasts,HitLeft(:,(Deflections==0)),VWeights,[0 .8 0],'filled')
set(gca,'FontSize',15)
legend(gca,'off');
title('Mouse choices - Visual Only','FontSize',20)
ylabel('Percent choice','FontSize',20)
xlabel('<-LEFT      log(10*Contrast)      RIGHT->','FontSize',20)

%TAC
x = Deflections;
yR = HitRight((Contrasts==0),:);
yL = HitLeft((Contrasts==0),:);
startPointsR = [ 0.15 , 0.35 , 60, 3];
startPointsL = [ 0.15 , 0.35 , -60, 3];
VWeights = TotalTrials((Contrasts==0),:);
ftacR = fit(x,yR',ft2,'Start', startPointsR, 'Weights', VWeights);
ftacL = fit(x,yL',ft2,'Start', startPointsL, 'Weights', VWeights);
TacCoeffsR = coeffvalues(fvisR);
TacCoeffsL = coeffvalues(fvisL);
%plot
figure();
hR = plot(ftacR,x,yR);
set(hR,'Color',[.8 0 0],'LineWidth',2)
hold on
hL = plot(ftacL,x,yL);
set(hL,'Color',[0 .8 0],'LineWidth',2)
ylim([0 1])
scatter(Deflections,HitRight((Contrasts==0),:),VWeights,[.8 0 0],'filled')
scatter(Deflections,HitLeft((Contrasts==0),:),VWeights,[0 .8 0],'filled')
set(gca,'FontSize',15)
legend(gca,'off');
title('Mouse choices - Tactile Only','FontSize',20)
ylabel('Percent choice','FontSize',20)
xlabel('<-LEFT      Deflection Intensity      RIGHT->','FontSize',20)