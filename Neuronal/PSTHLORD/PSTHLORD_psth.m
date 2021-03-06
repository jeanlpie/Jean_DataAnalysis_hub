%% Make raster
figure()
rast=[];
rastcon = [];
raw_psth={};
subplot(2,2,[1,3])
hold on %Allow multiple plots on the same graph
colorvec = colormap(winter(size(STS,2))); %colormap
tempsms = win(1):win(2);
frpertrial = nan(length(EV_TS),length(tempsms));
Gauss_width = max([11 6*FilterSize+1]);
kernel      = normpdf(-floor(Gauss_width/2):floor(Gauss_width/2),0,FilterSize);

for ch=1:size(STS,2) %for each channel 
    for k=1:length(EV_TS) %number of trials
    app = find(STS{1,ch}>=(EV_TS(k)+win(1)*1E3) & STS{1,ch}<=(EV_TS(k)+win(2)*1E3)); %find spikes ts in window
    raw_psth{k} = STS{1,ch}(app)-EV_TS(k); %make times relative to event
    t = raw_psth{k}; 
    rast = [rast, t]; %Stores ALL relative spiketimes for PSTH
    %make it in timebase 1ms then convolve w gaussian.
    frpertrial(k,:) = conv( histc(raw_psth{k}/1E3, tempsms) , kernel, 'same' ) / (size(STS,2) * 1E-3); %can use valid for only valid part of conv
    
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
Zpsth = ( psth - mu_bsl ) ./ sigma_bsl;
% filename=strcat('z-score_', condition_name, '.mat')
% save(filename, 'Zpsth', 'psth')

%Average the trial-by-trial firing rate
FRmean = mean(frpertrial);
FRstd = std(frpertrial); 
FRsem = FRstd/sqrt(length(EV_TS));
subplot(2,2,4)
shadedErrorBar(tempsms,FRmean,FRsem,'lineProps',{'Color',[.7 0 .7],'LineWidth',2},'patchSaturation',0.1)
xlim(win)
ylabel('Firing Rate (Hz)')

%Old SDF (smoothing psth)
% A=[];
% A(:,1) = (edges+BinSize*1e3/2)*1e-6;
% if boolZscored
%     A(:,2) = Zpsth;
% else
%     A(:,2) = psth;
% end
% [B] = msdf(A,'Gauss',0.5);
% subplot(2,2,4)
% plot(B(:,1),B(:,2))
% xlabel('Time (sec)') %Label x-axis
% if boolZscored
%     ylabel('Z-score') %Label y-axis
% else
%     ylabel('Hz') %Label y-axis
% end
% xlim([win(1) win(2)]*1e-3);

% filename=strcat('z-score', condition_name, '.pdf');
% fig = gcf;
% fig.PaperPositionMode = 'auto'
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(gcf,filename,'-dpdf')
% pause(0.1)
% close