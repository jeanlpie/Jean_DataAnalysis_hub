%% PSTH and Z-SCORE OF BASELINE WITHOUT PLOT
rast=[];
raw_psth={};
hold on %Allow multiple plots on the same graph
colorvec = colormap(winter(size(STS,2))); %colormap

for ch=1:size(STS,2) %for each channel 
    for k=1:length(EV_TS) %number of trials
    app = find(STS{1,ch}>=(EV_TS(k)+win(1)*1E3) & STS{1,ch}<=(EV_TS(k)+win(2)*1E3)); %find spikes ts in window
    raw_psth{k} = STS{1,ch}(app)-EV_TS(k); %make times relative to event
    t = raw_psth{k};
    rast = [rast, t]; %Stores ALL relative spiketimes for PSTH
    end
end


edges=(win(1):BinSize:win(2))*1e3;
n=histc(rast,edges);
psth=n./(length(EV_TS)*size(STS,2)*BinSize*1e-3); %normalize for bin duration, number of events and number of channels (since it's a grand total psth)

mu_bsl = mean(psth);
sigma_bsl = std(psth); 