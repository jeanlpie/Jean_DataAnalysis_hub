%% Neuronal data analysis after UmbertoSort
%% Separate data into clusters (Don't run if already have clustered spikes_ts and waveforms)
clu = 0;
for ch=1:size(clustersToKeep,2) 
   vec = clustersToKeep{1,ch};
   for k=1:length(vec)
       clu = clu+1;
       spikes_clu_ts{1,clu}         = spikes_ts{1,ch}(spikes_clusters{1,ch}==vec(k));          
       spikes_clu_waveforms{1,clu}  = spikes_waveforms{1,ch}((spikes_clusters{1,ch}==vec(k)),:);
       clus.channel(clu)            = ch;
   end
end
clear clu vec ch k

%% Get events (Don't run if you already have events)
%CHANGE THIS :
FileName='D:\cheetahData\J011_recs\2018-05-19_17-58-01_J011_1st\Events.nev'; %!! CAREFUL!! 
%%%%%
FieldSelectionFlags = [1 0 1 0 0];
[events_ts, events_ttl] = Nlx2MatEV( FileName, FieldSelectionFlags, 0 , 1);
clear FieldSelectionFlags FileName

%% Create events_str (Don't run if already have events_str)
%CAREFUL ! Use the proper TTL to event-type conversion.
events_str=cell(1,length(events_ttl));
%
events_str(events_ttl==1)= {'blank'};
events_str(events_ttl==2)= {'blank_opto'};
events_str(events_ttl==3)= {'visual_L'};
events_str(events_ttl==4)= {'visual_L_opto'};
events_str(events_ttl==5)= {'tactile_L'};
events_str(events_ttl==6)= {'tactile_L_opto'};
events_str(events_ttl==7)= {'multi_L'};
events_str(events_ttl==8)= {'multi_L_opto'};
events_str(events_ttl==9)= {'catch'};
events_str(events_ttl==10)= {'catch_opto'};
events_str(events_ttl==11)= {'visual_R'};
events_str(events_ttl==12)= {'visual_R_opto'};
events_str(events_ttl==13)= {'tactile_R'};
events_str(events_ttl==15)= {'tactile_R_opto'}; %%there was a mistake here
events_str(events_ttl==16)= {'multi_R'};    %here too
events_str(events_ttl==17)= {'multi_R_opto'}; %%till here
%
events_str(events_ttl==23)= {'LFRend'};
events_str(events_ttl==25)= {'checker'};
events_str(events_ttl==31)= {'endtrial'};
%
events_str(events_ttl==32768)= {'lickRight'};
events_str(events_ttl==16384)= {'lickLeft'};




%% Split or sub-select from clusters
STS     = spikes_clu_ts;
EV_TS   = events_ts;
EV_STR  = events_str;
SWF     = spikes_clu_waveforms;

% STS     = spikes_ts;
% EV_TS   = events_ts;
% EV_STR  = events_str;
% SWF     = spikes_waveforms;

%% Cluster quality assessment
%Select which clusters to analyse here:
CluSel = [2];%<--HERE
%%%
for clu = CluSel
    figure()
    % Waveform + average plots 
    subplot(2,2,[1 3])
    plot(SWF{1,clu}','Color',[.5 .5 .5])
    hold on
    plot(mean(SWF{1,clu}),'Color',[.8 .1 .1],'LineWidth',3)
    xlim([0 64])
    ylim([-600 900])
    title('Waveform');ylabel('uV');xlabel('sample number')
    % Inter-Spike-Interval histogram
    subplot(2,2,2)
    ISI=diff(STS{1,clu}); ISI=ISI/1000; ISI=ISI(ISI<350);
    histogram(ISI,'BinWidth',1)
    xlim([0 350])
    xlabel('ms')
    title('ISI')
    % Inter-Spike-Interval histogram
    subplot(2,2,4)
    ISI=diff(STS{1,clu}); ISI=ISI/1000; ISI=ISI(ISI<20);
    histogram(ISI,'BinWidth',1)
    xlim([0 20])
    xlabel('ms')
    title('ISI')
    %Title stating cluster number
    suptitle(strcat('Cluster number:',num2str(clu)))
end

%% Raster-plots and psths
%Now it's about aligning neural activity to an event. Select which event to
%align to :
    e1 = 'visual'; %Event to align to

% %
% for clu =  %12%CluSel
%     [~,~]=make_psth5(EV_TS,STS(clu),EV_STR, -1400000,2000000,50000,e1,1);
%     %Times in microseconds : pre, post, and bin.
%     suptitle(strcat('Cluster n°',num2str(clu),'-',e1));
% end
%Uncomment if you'd like to plot all selected clusters together
[~,~]=make_psth5(EV_TS,STS(CluSel),EV_STR, -500000,1600000,50000,e1,0);
%title(e1)
%ylim([0 14])
%save 24 16 8 2

%% Scatter-Response
%Now it's about comparing the neural response of two different events.
%Select which ones:
    e1 = 'visual_R2';
    e2 = 'visual_R1';
%
numclu = length(STS);
Resp=nan(numclu,2);
AR=nan(numclu,2);
Sil=[];
Labl=[];
BaseL = [];
EventResp = [];
    
% PSTH per cluster. For tactile and optotactile.
for clu = 32:82;%numclu %CluSel %1:numclu
    % for e1
    [psth_t_pre, ~] = mpsth(STS{1,clu}./1000000,EV_TS(strcmp(EV_STR,e1))./1000000,'pre',1505,'post',-5);
    [psth_t_post, ~] = mpsth(STS{1,clu}./1000000,EV_TS(strcmp(EV_STR,e1))./1000000,'pre',1,'post',501);  
    NumT = sum(strcmp(EV_STR,e1));
    %AR(clu,1) = Resp(clu,1) / ( (sum(psth_t_post(:,2)) / (NumT*0.5))  +  (sum(psth_t_pre(:,2)) / (NumT*1.5)) );                                %firing rate post.
    BaseL(clu,1) = (sum(psth_t_pre(:,2)) / (NumT*1.5));
    EventResp(clu,1) = (sum(psth_t_post(:,2)) / (NumT*0.5));

    % for tactile_opto
    [psth_to_pre, ~] = mpsth(STS{1,clu}./1000000,EV_TS(strcmp(EV_STR,e2))./1000000,'pre',1505,'post',-5);
    [psth_to_post, ~] = mpsth(STS{1,clu}./1000000,EV_TS(strcmp(EV_STR,e2))./1000000,'pre',1,'post',501); 
    NumTO = sum(strcmp(EV_STR,e2));
    %AR(clu,2)= Resp(clu,2) / ( (sum(psth_to_post(:,2)) / (NumTO*0.5))  +  (sum(psth_to_pre(:,2)) / (NumTO*1.5)) );
    BaseL(clu,2) = (sum(psth_to_pre(:,2)) / (NumTO*1.5));
    EventResp(clu,2) = (sum(psth_to_post(:,2)) / (NumTO*0.5));


    Labl{1,clu}=num2str(clu);
    Sil(clu) = ( (sum(psth_to_post(:,2)) / (NumTO*0.5)) - (sum(psth_t_post(:,2)) / (NumT*0.5)) ) / ( (sum(psth_to_post(:,2)) / (NumTO*0.5)) + (sum(psth_t_post(:,2)) / (NumT*0.5)) );      
end
Resp = EventResp - BaseL;

% Plot
figure()
scatter(Resp(:,1),Resp(:,2))
text(Resp(:,1),Resp(:,2),Labl,'FontSize',10,'Color',[.5 0 0])
title('Response-Baseline for each cluster')
ylabel(strcat(e2,'(Hz)'))
xlabel(strcat(e1,'(Hz)'))
line([min(min(Resp)) max(max(Resp))],[min(min(Resp)) max(max(Resp))])
    
%Plot opto APr versus no-opto APr
%     figure()
%     scatter(AR(:,1),AR(:,2))
%     title('(R-B)/(R+B) for each cluster')
%     ylabel(strcat(e2,'(Hz)'))
%     xlabel(strcat(e1,'(Hz)'))
%     line([min(min(AR)) max(max(AR))],[min(min(AR)) max(max(AR))])

%     figure()
%     scatter(clus.lpos,Sil)

figure()
boxplot(Sil,clus.lpos)

%% Cross-correlation
%Compare the following clusters:
    clus1 = 13;
    clus2 = 14;
%
figure()
%
subplot(1,2,1)
[tsOffsets] = crosscorrelogram(STS{clus1}, STS{clus2}, [-40e3 40e3]);
histogram(tsOffsets,'BinWidth',1e3)
title(strcat('Crosscorr: ',num2str(clus2),' rel to ',num2str(clus1)))
%
subplot(1,2,2)
plot(mean(SWF{clus1}))
hold on
plot(mean(SWF{clus2}))

%% Classify as putative IN or PyrN
fs = 32552;
pk2thr = nan(1,length(SWF));
for clu = 1:length(SWF)
meanwv = mean(SWF{1,clu},1);
[~, t_max] = max(meanwv);
[~, t_min] = min(meanwv(t_max+1:end));
pk2thr(clu) = t_min/fs*1000;
meanwv=[];
t_max = [];
t_min = [];
end
clear meanwv t_max t_min

histogram(pk2thr,'BinWidth',0.01)

% Classify the clusters as FSint if <=0.35 and pPyrN if >=0.45ms
% Unclassified : 0
% FSInt : 1
% pPyrN : 2
prompt = {'Enter lower cutoff:','Enter higher cutoff:'};
dlg_title = 'Select threshold for classification';
num_lines = 1;
defaultans = {'0.35','0.45'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

clu_class = zeros(1,length(SWF));
clu_class(pk2thr<=0.35) = 1; 
clu_class(pk2thr>=0.45) = 2; 

spikes_clu_ts_fsi = STS(clu_class==1);
spikes_clu_ts_ppn = STS(clu_class==2);

spikes_clu_waveforms_fsi = SWF(clu_class==1); 
spikes_clu_waveforms_ppn = SWF(clu_class==2);

%Normalize waveforms means
for k=1:length(spikes_clu_waveforms_fsi)
    A = mean(spikes_clu_waveforms_fsi{1,k});
    A = A./max(A);
    WFMN_FSI(:,k)=A;
end
for k=1:length(spikes_clu_waveforms_ppn)
    A = mean(spikes_clu_waveforms_ppn{1,k});
    A = A./max(A);
    WFMN_PPN(:,k)=A;
end

%plot waveforms
figure()
plot(WFMN_FSI,'r')
hold on
plot(WFMN_PPN,'b')
% for k=1:length(spikes_clu_waveforms_ppn)
%    plot(mean(spikes_clu_waveforms_ppn{1,k}),'b') 
%    hold on
% end
% for k=1:length(spikes_clu_waveforms_fsi)
%    plot(mean(spikes_clu_waveforms_fsi{1,k}),'r') 
%    hold on
% end
