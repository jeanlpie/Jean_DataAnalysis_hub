%% Raster, PSTH and SDF 
%This code allows the user to select the desired data to analyze, splits it
%according to given conditions and plots rasters per cluster and
%grand-total psth and sdf.
%The desired parameters must be all specified first, then the code runs.
%You need to load :
%1- Neuronal data
%2- Behavioral data
%3- Events data
%% Select original data 
%Select the timestamp data to use (cell). Cells separate channels/clusters.
STS     = spikes_clu_ts; 
%Select the events timestamps (vector).
EV_TS   = events_ts; 
%Select the corresponding event names (cell with strings).
EV_STR  = events_str;

%% Subselect data
%Direct choice:
booluseDirect = true;
CluSelDirect = [1] %2 6 9 34[5 9 12 13 18 19 21];%[2 6 12 21 24 27];

%Probe area
boolAllAreas = false; %Use all available areas?
V1     = 0;
PPC    = 1;
RL     = 0;
Cg1    = 0;

%Cluster depth. Subselects only clusters which depths are in interval.
boolAllDepths     = true; %Use all depths ?
SS.Depth_interval = [0 1200]; %inclusive

%% Select conditions to split data
%Subselect or skip (if no trialData) ?
boolFilt = true;

%Trial-wise:
%TrialType (Modality)
trialType_catch     = 0;
trialType_tactile   = 0;
trialType_visual    = 1;
trialType_multi     = 0;

%Condition-number (For correct right V4, T9, M13)
SelConds            = 9; %1:32 %32 possible conditions [17 18 20 21]

%Stimulus Side (L vs R)
%ONLY USE FOR CONGRUENT SESSIONS (for now)
stimSide_L          = 0;
stimSide_R          = 1;
stimSide_N          = 0;

%Mouse choice resp.
firstRespLFR_L      = 0; 
firstRespLFR_R      = 1; 
firstRespLFR_none   = 0;

%Trial Outcome
correctResponse     = 1;
firstIncorrect      = 0;
noResponse          = 0;

%Optogenetics
laserON             = 0;
laserOFF            = 1;

%Other event sub-select (not stimuli):
LFR_End             = 0;
lickLeft            = 0;
lickRight           = 0;
%Epoch sub-select:
only_StimEpoch      = 0;
only_LFR            = 0;
only_ITI            = 0;


% Compose name: AREA_Modality_Side_Condition_Outcome_Opto
condition_name = strcat(SS.Area.names{SS.Area.bools}, '_');
if booluseDirect
    allOneString = sprintf('%.0f-' , CluSelDirect);
    allOneString = allOneString(1:end-1);
    condition_name = strcat(condition_name, allOneString);
    condition_name = strcat(condition_name, '_');
end
if trialType_catch; condition_name = strcat(condition_name, 'C'); end
if trialType_tactile; condition_name = strcat(condition_name, 'T'); end
if trialType_visual; condition_name = strcat(condition_name, 'V'); end
if trialType_multi; condition_name = strcat(condition_name, 'M'); end
condition_name = strcat(condition_name, '_');
if stimSide_L; condition_name = strcat(condition_name, 'L'); end
if stimSide_R; condition_name = strcat(condition_name, 'R'); end
if stimSide_N; condition_name = strcat(condition_name, 'N'); end
allOneString = sprintf('%.0f-' , SelConds);
allOneString = allOneString(1:end-1);
condition_name = strcat(condition_name, allOneString, '_');

if correctResponse; condition_name = strcat(condition_name, 'Co'); end
if firstIncorrect; condition_name = strcat(condition_name, 'In'); end
if noResponse; condition_name = strcat(condition_name, 'Nr'); end
condition_name = strcat(condition_name, '_');

if laserON; condition_name = strcat(condition_name, 'On'); end
if laserOFF; condition_name = strcat(condition_name, 'Off'); end

disp(condition_name)

%% Select analysis parameters
%Peri-Event Window for all plots
win             = [ -800 , 1200 ]; %in ms

%Bin Size for PSTH
BinSize         = 50; %in ms

%Baseline normalized for SDF
boolBslnNorm    = 0;

%Z-scored for SDF
boolZscored     = 0;

%Gaussian filter for SDF
boolGaussFilt   = 1;
FilterSize      = 20; %in ms

%% Select graphical plotting parameters
%Make new figure
boolMakeNewFigure   = 1;
%Color options 
Rasterkleurmap      = winter;
%If sorted by trialType
kleur{1} = [[0.4 0.4 0.4];[0.4 0.7 1];[1 0.7 0.5];[0.7 0.2 0.7]];%CTVM
kleur{2} = [[0 1 0.3];[1 0.3 0];[0.55 0.55 0.55]];%LRN (stimside)
kleur{3} = [[0 0.7 0];[0.7 0 0];[0.6 0.6 0.6]];%LRN (resp)
kleur{4} = [[0.85 0.8 0];[0 0 0];[0.5 0.5 0.5]];%outcome(C/I/NR)
kleur{4} = [[0 0.2 1];[0.7 0.7 0.9]];%Opto(ON/OFF)

%Sort by
%depth

%% Actual code
%Subselect data
CluSel = 1:length(STS); CluSel4 = 1:length(STS); CluSel2 = CluSel; CluSel1 = CluSel;
if boolAllAreas
    SS.Area.names = 'all';
else
SS.Area.names = {'V1','PPC','RL','Cg1'}; SS.Area.bools = logical([V1 PPC RL Cg1]); 
CluSel1 = CluSel1(ismember(clus.area,SS.Area.names(SS.Area.bools)));
end
if boolAllDepths
    SS.Depth_interval = [-inf +inf];
else
    CluSel2 = find(clus.zpos<=SS.Depth_interval(2) & clus.zpos>=SS.Depth_interval(1));
end
CluSel4 = CluSel1(ismember(CluSel1,CluSel2));
if booluseDirect
    CluSel = CluSelDirect(ismember(CluSelDirect,CluSel4));
else
    CluSel = CluSel4;
end
STS = STS(CluSel);
clear V1 PPC RL Cg1 CluSel* boolAllDepths boolAllAreas

%% clear selected events
idx_stims = find(ismember(events_ttl,[1:32 [1:32]+16384 [1:32]+32768])==1); %indexes of stims. can do with str too
%Subselect trials
if boolFilt
%modality
Filt1 = zeros(1,length(EV_TS));
if trialType_catch
  SCAT=zeros(1,length(EV_TS));
  SCAT(idx_stims(trialData.trialType=='C'))=1;
  Filt1 = Filt1|SCAT;
  %Filt1 = Filt1|strncmp(events_str,'catch',5)  ;
end
if trialType_tactile
  STAC=zeros(1,length(EV_TS));
  STAC(idx_stims(trialData.trialType=='T'))=1;
  Filt1 = Filt1|STAC;
  %Filt1 = Filt1|strncmp(events_str,'tactile',6) ;
end
if trialType_visual
  SVIS=zeros(1,length(EV_TS));
  SVIS(idx_stims(trialData.trialType=='V'))=1;
  Filt1 = Filt1|SVIS;
%   Filt1 = Filt1|strncmp(events_str,'visual',6) ;
end
if trialType_multi
  SMUL=zeros(1,length(EV_TS));
  SMUL(idx_stims(trialData.trialType=='M'))=1;
  Filt1 = Filt1|SMUL;
  %Filt1 = Filt1|strncmp(events_str,'multi',5) ;
end
%stimulus side
Filt2 = zeros(1,length(EV_TS));
if ~isfield(trialData,'stimSide') %for the rec2 version
    trialData.stimSide = trialData.stimTACSide;
    %only for congruency, otherwise visual overwrites.
    trialData.stimSide(strcmp(trialData.stimVISSide,'R'))={'R'};
    trialData.stimSide(strcmp(trialData.stimVISSide,'L'))={'L'};
end
    if stimSide_L
        SL = zeros(1,length(EV_TS));
        SL(idx_stims(strcmp(trialData.stimSide,'L')))=1; %as the rank should be the same.
        Filt2 = Filt2|SL;
    end
    if stimSide_R
        SR = zeros(1,length(EV_TS));
        SR(idx_stims(strcmp(trialData.stimSide,'R')))=1;
        Filt2 = Filt2|SR;
    end
    if stimSide_L
        SN = zeros(1,length(EV_TS));
        SN(idx_stims(strcmp(trialData.stimSide,'N')))=1;
        Filt2 = Filt2|SN;
    end
%Mouse response
Filt3 = zeros(1,length(EV_TS));
if firstRespLFR_L
    RL = zeros(1,length(EV_TS));
    RL(idx_stims(trialData.firstRespLFR=='L'))=1;
    Filt3 = Filt3|RL;
end
if firstRespLFR_R
    RR = zeros(1,length(EV_TS));
    RR(idx_stims(trialData.firstRespLFR=='R'))=1;
    Filt3 = Filt3|RR;
end
if firstRespLFR_none
    RN = zeros(1,length(EV_TS));
    RN(idx_stims(isnan(trialData.firstRespLFR)))=1;
    Filt3 = Filt3|RN;
end
%Trial Outcome
Filt4 = zeros(1,length(EV_TS));
if correctResponse
    CR = zeros(1,length(EV_TS));
    CR(idx_stims(trialData.correctResponse==1))=1;
    Filt4 = Filt4|CR;
end
if firstIncorrect    
    FI = zeros(1,length(EV_TS));
    FI(idx_stims(trialData.firstIncorrect==1))=1;
    Filt4 = Filt4|FI;
end
if noResponse        
    NR = zeros(1,length(EV_TS));
    NR(idx_stims(trialData.noResponse==1))=1;
    Filt4 = Filt4|NR;
end
%Optogenetics
Filt5 = zeros(1,length(EV_TS));
if laserON     
    LON = zeros(1,length(EV_TS));
    LON(idx_stims(trialData.laserON==1))=1;
    Filt5 = Filt5|LON;    
end
if laserOFF       
    LOFF = zeros(1,length(EV_TS));
    LOFF(idx_stims(trialData.laserON==0))=1;
    Filt5 = Filt5|LOFF; 
end
%Condition-number
Filt6 = zeros(1,length(EV_TS));
    FCON = zeros(1,length(EV_TS));
    FCON(idx_stims(ismember(trialData.conditionnumber,SelConds)))=1;
    Filt6 = Filt6|FCON;  
%Merge filters
Filt = Filt1&Filt2&Filt3&Filt4&Filt5&Filt6;
STIMS_TS = EV_TS(Filt);
%Subselect other event or just stims ?(
if ~(LFR_End||lickLeft||lickRight) 
    %By default just use stims
    EV_TS = STIMS_TS;
else
    %Use other events like licks or LFRend in the subselected trials.
    if lickLeft
    LL_TS = EV_TS(find(strcmp(EV_STR,'lickLeft')));
    %Only for trying first lick in bout
    DTS = diff([0 LL_TS]);
    LL_TS = LL_TS(DTS>300000);
    %
    EV_TS = LL_TS;
    end
    if lickRight
    LR_TS = EV_TS(find(strcmp(EV_STR,'lickRight')));
    %Only for trying first lick in bout
    DTS = diff([0 LR_TS]);
    LR_TS = LR_TS(DTS>300000);
    %
    EV_TS = LR_TS;
    end
end

%Subselect epoch

if only_ITI
    ITIEV = [];
    for k=1:length(STIMS_TS)
        ITIEV = [ITIEV EV_TS(EV_TS>(STIMS_TS(k)-2.5*1E6) & EV_TS<(STIMS_TS(k)-10))]; 
    end
    EV_TS = ITIEV;
end

%Subselect

end



%Choose color

%% Make raster
figure()
rast=[];
raw_psth={};
subplot(2,2,[1,3])
hold on %Allow multiple plots on the same graph
colorvec = colormap(winter(size(STS,2))); %colormap

for ch=1:size(STS,2) %for each channel 
    for k=1:length(EV_TS) %number of trials
    app = find(STS{1,ch}>=(EV_TS(k)+win(1)*1E3) & STS{1,ch}<=(EV_TS(k)+win(2)*1E3)); %find spikes ts in window
    raw_psth{k} = STS{1,ch}(app)-EV_TS(k); %make times relative to event
    t = raw_psth{k};
    rast = [rast, t]; %Stores ALL relative spiketimes for PSTH
    
    for i = 1:length(t) %Loop through each spike time
        line([t(i) t(i)]*1e-6, [k+((ch-1)*length(EV_TS))-1 k+((ch-1)*length(EV_TS))], 'Color', colorvec(ch,:),'LineWidth',2) %Create a tick mark at x = t1(i) with a height of 1
    end
    end
    line ([win(1)*1E-3 ((win(2)-win(1))/6+win(1))*1E-3],[ch*length(EV_TS) ch*length(EV_TS)],'Color', colorvec(ch,:),'LineWidth', 1); %Line separates channels
end
xlim([win(1) win(2)]*1e-3);
%Make PSTH
subplot(2,2,2);
edges=(win(1):BinSize:win(2))*1e3;
n=histc(rast,edges);
psth=n./(length(EV_TS)*size(STS,2)*BinSize*1e-3); %normalize for bin duration, number of events and number of channels (since it's a grand total psth)
bar(edges*1e-6,psth);
xlabel('Time (sec)') %Label x-axis
ylabel('Hz') %Label x-axis
xlim([win(1) win(2)]*1e-3);

%Get Z-scores
Zpsth=zscore(psth)
filename=strcat('z-score_', condition_name, '.mat')
save(filename, 'Zpsth', 'psth')

%Make SDF
A=[];
A(:,1) = (edges+BinSize*1e3/2)*1e-6;
A(:,2) = psth;
[B] = msdf(A,'Gauss',2);
subplot(2,2,4)
plot(B(:,1),B(:,2))
xlabel('Time (sec)') %Label x-axis
ylabel('Hz') %Label y-axis
xlim([win(1) win(2)]*1e-3);

filename=strcat('z-score', condition_name, '.pdf');
fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf,filename,'-dpdf')
pause(0.1)
close

%% Plot several firing frequencies in one plot

fc = load('z-score_PPC_28_T_R9_Co_Off')
red_line = fc.Zpsth

fc = load('z-score_PPC_28_V_R4_Co_Off')
blue_line = fc.Zpsth

fc = load('z-score_PPC_28_M_R13_Co_Off.mat')
purple_line = fc.Zpsth

labels = {'TR', 'VR', 'MR'} % Red blue purple labels
filename=strcat('frequencies_newplot.pdf');

%Make SDF with Z-scores
figure()
[B] = msdf([((edges+BinSize*1e3/2)*1e-6)', red_line'],'Gauss',2);
plot(B(:,1),B(:,2), 'r')
hold on
[B] = msdf([((edges+BinSize*1e3/2)*1e-6)', blue_line'],'Gauss',2);
plot(B(:,1),B(:,2), 'b')
[B] = msdf([((edges+BinSize*1e3/2)*1e-6)', purple_line'],'Gauss',2);
plot(B(:,1),B(:,2), 'color', [0.49 0.18 0.45])

% title('Correct')
legend(labels, 'Location', 'NorthWest')
xlabel('Time (sec)') %Label x-axis
ylabel('psth z-score') %Label y-axis
set(findall(gcf,'-property','FontSize'),'FontSize',18)
xlim([win(1) win(2)]*1e-3);
ylim([-1 1]*1.5);

fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(gcf,filename,'-dpdf')
pause(0.1)



% boolBslnNorm    = 0;
% boolZscored     = 0;
% boolGaussFilt   = 1;
% FilterSize      = 25; %in ms
