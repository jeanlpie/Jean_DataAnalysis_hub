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