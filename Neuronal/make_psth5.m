function [raw_psth,psth]=make_psth5(events_ts,spikes_ts,events_string, t_pre,t_post,psth_bin,event_str,startswith)
% psth all trials and all channels. 
% All channels per trial are combined into the same line in the raster.
% Channels are color coded from green to blue. To deactivate just erase the
% color part in line 24.
% This also sorts by channel
% And you can choose which event.

edges=t_pre:psth_bin:t_post; %in micro-seconds
lines = 1; rows = 2; subplot1 = 1; subplot2 = 2; %default
figure()

rast=[];
raw_psth={};
subplot(lines,rows,subplot1);
hold on %Allow multiple plots on the same graph
colorvec = colormap(winter(size(spikes_ts,2))); %colormap

if startswith
   events_ts = events_ts( strncmpi(events_string,event_str,length(event_str)) );
else
events_ts = events_ts(strcmp(events_string,event_str));
end

%ERASE THIS SHIT NOW:
%events_ts = events_ts(1+0*18:18+0*18);
%END

for ch=1:size(spikes_ts,2) %for each channel 
  
    for k=1:length(events_ts) %number of trials
    
    %find spikes_ts inside the window around each event. 
   
    app = find(spikes_ts{1,ch}>=(events_ts(k)+t_pre) & spikes_ts{1,ch}<=(events_ts(k)+t_post)); %find spikes ts in window
    raw_psth{k} = spikes_ts{1,ch}(app)-events_ts(k); %make times relative to event
    t = raw_psth{k};
    rast = [rast, t]; %Stores ALL relative spiketimes for PSTH
    for i = 1:length(t) %Loop through each spike time
        line([t(i) t(i)]*1e-6, [k+((ch-1)*length(events_ts))-1 k+((ch-1)*length(events_ts))], 'Color', colorvec(ch,:),'LineWidth',2) %Create a tick mark at x = t1(i) with a height of 1
    end
    end
    line ([t_pre*1e-6 ((t_post-t_pre)/6+t_pre)*1e-6],[ch*length(events_ts) ch*length(events_ts)],'Color', colorvec(ch,:),'LineWidth', 1); %Line separates channels
end
xlim([t_pre t_post]*1e-6);
n=histc(rast,edges);
psth=n./(length(events_ts)*size(spikes_ts,2)*psth_bin*1e-6); %normalize for bin duration, number of events and number of channels (since it's a grand total psth)
subplot(lines,rows,subplot2);
bar(edges*1e-6,psth);
xlabel('Time (sec)') %Label x-axis
ylabel('Hz') %Label x-axis
xlim([t_pre t_post]*1e-6);
