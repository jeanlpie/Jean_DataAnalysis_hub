%% Select conditions to split data

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
CluSelDirect = [13] %2 6 9 34[5 9 12 13 18 19 21];%[2 6 12 21 24 27];

%Probe area
boolAllAreas = false; %Use all available areas?
V1     = 0;
PPC    = 1;
RL     = 0;
Cg1    = 0;

%Cluster depth. Subselects only clusters which depths are in interval.
boolAllDepths     = true; %Use all depths ?
SS.Depth_interval = [0 1200]; %inclusive

%Subselect or skip (if no trialData) ?
boolFilt = true;

%Trial-wise:
%TrialType (Modality)
trialType_catch     = 1;
trialType_tactile   = 1;
trialType_visual    = 1;
trialType_multi     = 1;

%Condition-number (For correct right V4, T9, M13)
SelConds            = 1:32; %1:32 %32 possible conditions [17 18 20 21]

%Stimulus Side (L vs R)
%ONLY USE FOR CONGRUENT SESSIONS (for now)
stimSide_L          = 1;
stimSide_R          = 1;
stimSide_N          = 1;

%Mouse choice resp.
firstRespLFR_L      = 1; 
firstRespLFR_R      = 1; 
firstRespLFR_none   = 1;

%Trial Outcome
correctResponse     = 1;
firstIncorrect      = 1;
noResponse          = 1;

%Optogenetics
laserON             = 1;
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
% condition_name = strcat(SS.Area.names{SS.Area.bools}, '_');
% if booluseDirect
%     allOneString = sprintf('%.0f-' , CluSelDirect);
%     allOneString = allOneString(1:end-1);
%     condition_name = strcat(condition_name, allOneString);
%     condition_name = strcat(condition_name, '_');
% end
% if trialType_catch; condition_name = strcat(condition_name, 'C'); end
% if trialType_tactile; condition_name = strcat(condition_name, 'T'); end
% if trialType_visual; condition_name = strcat(condition_name, 'V'); end
% if trialType_multi; condition_name = strcat(condition_name, 'M'); end
% condition_name = strcat(condition_name, '_');
% if stimSide_L; condition_name = strcat(condition_name, 'L'); end
% if stimSide_R; condition_name = strcat(condition_name, 'R'); end
% if stimSide_N; condition_name = strcat(condition_name, 'N'); end
% allOneString = sprintf('%.0f-' , SelConds);
% allOneString = allOneString(1:end-1);
% condition_name = strcat(condition_name, allOneString, '_');
% 
% if correctResponse; condition_name = strcat(condition_name, 'Co'); end
% if firstIncorrect; condition_name = strcat(condition_name, 'In'); end
% if noResponse; condition_name = strcat(condition_name, 'Nr'); end
% condition_name = strcat(condition_name, '_');
% 
% if laserON; condition_name = strcat(condition_name, 'On'); end
% if laserOFF; condition_name = strcat(condition_name, 'Off'); end
% 
% disp(condition_name)

%% Select analysis parameters
%Peri-Event Window for all plots
win             = [ -2000 , 0 ]; %in ms

%Bin Size for PSTH
BinSize         = 50; %in ms

%Baseline normalized for SDF
boolBslnNorm    = 0;

%Z-scored for SDF
boolZscored     = 0;

%Gaussian filter for SDF
boolGaussFilt   = 1;
FilterSize      = 50; %in ms

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
