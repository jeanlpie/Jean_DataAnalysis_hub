%% Checks final task behavior with full-info plot
%% Cut trials
numtrials = length(trialData.trialEnd)-3;
tfields = fieldnames(trialData);
for k=1:length(tfields) %clean trials that didnt happen
        if (isnumeric( trialData.(tfields{k}) ) || iscell( trialData.(tfields{k}) )) && length(trialData.(tfields{k}))>numtrials
            trialData.(tfields{k}) = trialData.(tfields{k})(1:numtrials);
        end
end
%Cut the end ?

%%
%Correct if unimodal trials still have a value for the intensity of the non-shown modality.
%Actually pie_task only saves what will be shown. Zero values are therefore nans.
%Convert nans into zeros:
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;

% %Contrasts in log ?
% trialData.stimContrast = log10(trialData.stimContrast*1000);
% trialData.stimContrast(isinf(trialData.stimContrast)) = 0; %ATTENTION ! Contrast cannot be smaller or equal to 1.

%Make left negative
trialData.stimDeflection(trialData.leftCorrect==1) =  trialData.stimDeflection(trialData.leftCorrect==1)*-1;
trialData.stimContrast(trialData.leftCorrect==1)   =  trialData.stimContrast(trialData.leftCorrect==1)*-1;

%Supposing M,V and T trials had the same intensity values:
Deflections = unique(trialData.stimDeflection); Deflections(isnan(Deflections))=[];
Contrasts = unique(trialData.stimContrast); Contrasts(isnan(Contrasts))=[];

%Tweaky stuff
%  Contrasts([4 6])=[];
%  Deflections(4)=[];
% Deflections(6) = [];

NTac = length(Deflections);
NVis = length(Contrasts);

%  Compute percent right, left, and no response
for v=1:NVis
    for t=1:NTac
      TotalTrials(v,t) =  nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));  
      HitRight(v,t) = nansum(trialData.firstRespLFR=='R' & trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v)) / nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));
      HitLeft(v,t)  = nansum(trialData.firstRespLFR=='L' & trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v)) / nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));
      NoResp(v,t)   = nansum(isnan(trialData.firstRespLFR) & trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v)) / nansum(trialData.stimDeflection==Deflections(t) & trialData.stimContrast==Contrasts(v));
    end
end

%% Plot for VisOnly and TacOnly
figure()
hold on
plot(Contrasts,HitRight(:,(Deflections==0)),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Contrasts,HitLeft(:,(Deflections==0)),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Contrasts,NoResp(:,(Deflections==0)),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
ylim([0 1])
title('Mouse choices - Visual Only')
ylabel('Percent choice')
xlabel('<-LEFT      Contrast      RIGHT->')
legend({'LickRight','LickLeft','NoResponse'})

figure()
hold on
%For x-axis in degrees, replace Deflections by : atand(Deflections/311.68)
plot(Deflections,HitRight((Contrasts==0),:),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Deflections,HitLeft((Contrasts==0),:),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Deflections,NoResp((Contrasts==0),:),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
ylim([0 1])
title('Mouse choices - Tactile Only')
ylabel('Percent choice')
xlabel('<-LEFT      Deflection      RIGHT->')
legend({'LickRight','LickLeft','NoResponse'})

%% IMAGESCS
%figure('pos',[0 180 1400 500])
figure('Units','normalized','pos',[0.05 0.4 0.8 0.5])
colormap(summer)


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

%% New RT
figure()
%Visual
RT_V=cell(1,1);
for k = 1:length(trialData.stimStart)
    if trialData.trialType(k)=='V' %condition : visual
        RT_V{k} = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
        if ~isempty(RT_V{k})
            RT_V{k} = RT_V{k}(1);
            RT_V{k} = RT_V{k} - trialData.stimStart(k);
        end
    end
end
RT_V=cell2mat(RT_V);
%Tactile
RT_T=cell(1,1);
for k = 1:length(trialData.stimStart)
    if trialData.trialType(k)=='T' %condition: tactile
        RT_T{k} = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
        if ~isempty(RT_T{k})
            RT_T{k} = RT_T{k}(1);
            RT_T{k} = RT_T{k} - trialData.stimStart(k);
        end
    end
end
RT_T=cell2mat(RT_T);
%Multimodal
RT_M=cell(1,1);
for k = 1:length(trialData.stimStart)
    if trialData.trialType(k)=='M'
        RT_M{k} = trialData.lickTime{k}(trialData.lickTime{k}>=trialData.stimStart(k));
        if ~isempty(RT_M{k})
            RT_M{k} = RT_M{k}(1);
            RT_M{k} = RT_M{k} - trialData.stimStart(k);
        end
    end
end
RT_M=cell2mat(RT_M);
%Plot
subplot(3,1,1)
histogram(RT_V,'BinWidth',0.01)
title('Visual RT')
xlim([0 1])
subplot(3,1,2)
histogram(RT_T,'BinWidth',0.01)
title('Tactile RT')
xlim([0 1])
subplot(3,1,3)
histogram(RT_M,'BinWidth',0.01)
title('Multi RT')
xlim([0 1])
%% Sliding performance
size = 60;
numtrials = length(trialData.correctResponse);
perfect=(trialData.correctResponse==1 & trialData.trialType~='C')|(trialData.noResponse==1 & trialData.trialType=='C');%

total= ones(1,numtrials); %can insert extra condition here and in perfect.
perf = nan(numtrials,1);

for k=size:numtrials-1
        perf(k) = sum(perfect(k-size+1:k)) / sum(total(k-size+1:k));  %size
        hitsL(k) = sum(trialData.correctResponse(k-size+1:k)==1 & trialData.trialType(k-size+1:k)~='C' & trialData.leftCorrect(k-size+1:k)==1) / sum(trialData.trialType(k-size+1:k)~='C' & trialData.leftCorrect(k-size+1:k)==1);
        hitsR(k) = sum(trialData.correctResponse(k-size+1:k)==1 & trialData.trialType(k-size+1:k)~='C' & trialData.rightCorrect(k-size+1:k)==1) / sum(trialData.trialType(k-size+1:k)~='C' & trialData.rightCorrect(k-size+1:k)==1);
        fars(k) = sum(trialData.firstIncorrect(k-size+1:k)==1 & trialData.trialType(k-size+1:k)=='C') / sum(trialData.trialType(k-size+1:k)=='C');
        visperf(k)  = sum(trialData.correctResponse(k-size+1:k)==1 & trialData.trialType(k-size+1:k)=='V') / sum(trialData.trialType(k-size+1:k)=='V');
        tacperf(k)  = sum(trialData.correctResponse(k-size+1:k)==1 & trialData.trialType(k-size+1:k)=='T') / sum(trialData.trialType(k-size+1:k)=='T');
        multiperf(k)= sum(trialData.correctResponse(k-size+1:k)==1 & trialData.trialType(k-size+1:k)=='M') / sum(trialData.trialType(k-size+1:k)=='M');
end


figure()
subplot(2,2,1) %average perf
plot(perf,'LineWidth',3)
line([0 numtrials],[0.5 0.5],'Color','r')
ylim([0 1])
xlim([0 numtrials])
ylabel('Performance', 'FontSize', 28)
xlabel('Trial Number', 'FontSize', 28)
title(strcat('Average performance - ',sessionData.Mouse,' - ',sessionData.Date), 'FontSize', 28)

subplot(2,2,2) %Hits and FAs
plot(hitsL,'LineWidth',3,'Color',[0 .7 0])
hold on
plot(hitsR,'LineWidth',3,'Color',[.7 0 0])
plot(fars,'LineWidth',3,'Color',[.7 .7 .7])
line([0 numtrials],[0.5 0.5],'Color','r')
ylim([0 1])
xlim([0 numtrials])
legend({'HitLeft','HitRight','FARatio'},'Location','SouthWest')
xlabel('Trial Number', 'FontSize', 28)

subplot(2,2,3) %per modality
plot(visperf,'LineWidth',3,'Color',[1 0.7 0.5])
hold on
plot(tacperf,'LineWidth',3,'Color',[0.4 0.7 1])
plot(multiperf,'LineWidth',3,'Color',[0.7 0.2 0.7])
line([0 numtrials],[0.5 0.5],'Color','r')
ylim([0 1])
xlim([0 numtrials])
legend({'Visual','Tactile','Multi'},'Location','SouthWest')
xlabel('Trial Number', 'FontSize', 28)

%% Lick proba
licksC = [];
licksV = [];
licksT = [];
licksM = [];
weirdtrials = [];
for k=1:length(trialData.correctResponse)
    if ~isempty(trialData.lickTime{k,1})
            switch trialData.trialType(k)
                case 'C'
                    for i=1:length(trialData.lickTime{k,1})
                    lick = trialData.lickTime{k,1}(i) - trialData.stimStart(k);
                    licksC = [licksC lick];
%                         if lick>-0.1 && lick<0.1
%                             weirdtrials = [weirdtrials k];
%                         end
                    end
                case 'V'
                    for i=1:length(trialData.lickTime{k,1})
                    lick = trialData.lickTime{k,1}(i) - trialData.stimStart(k);
                    licksV = [licksV lick];
                    end
                case 'T'
                    for i=1:length(trialData.lickTime{k,1})
                    lick = trialData.lickTime{k,1}(i) - trialData.stimStart(k);
                    licksT = [licksT lick];
                        if lick>-0.1 && lick<0.1
                            weirdtrials = [weirdtrials k];
                        end
                    end
                case 'M'
                    for i=1:length(trialData.lickTime{k,1})
                    lick = trialData.lickTime{k,1}(i) - trialData.stimStart(k);
                    licksM = [licksM lick];
                    end                    
            end
    end
end         

figure()
edgiz = -2.5:0.1:1.2;

subplot(4,1,1)
hC = histc(licksC,edgiz); hC = hC/(sum(trialData.trialType=='C')*diff(edgiz(2:3)));
plot(edgiz,hC,'Color',[0.4 0.4 0.4],'LineWidth',1)
title('Lick rate in catch trials (licks per second per trial)')
xlim([-2.5 1])
ylim([0 6])

subplot(4,1,2)
hV = histc(licksV,edgiz); hV = hV/(sum(trialData.trialType=='V')*diff(edgiz(2:3)));
plot(edgiz,hV,'Color',[0.7 0.4 0],'LineWidth',1)
title('Lick rate in visual trials (licks per second per trial)')
xlim([-2.5 1])
ylim([0 6])

subplot(4,1,3)
hT = histc(licksT,edgiz); hT = hT/(sum(trialData.trialType=='T')*diff(edgiz(2:3)));
plot(edgiz,hT,'Color',[0 0.2 0.7],'LineWidth',1)
title('Lick rate in tactile trials (licks per second per trial)')
xlim([-2.5 1])
ylim([0 6])

subplot(4,1,4)
hM = histc(licksM,edgiz); hM = hM/(sum(trialData.trialType=='M')*diff(edgiz(2:3)));
plot(edgiz,hM,'Color',[0.7 0.2 0.7],'LineWidth',1)
title('Lick rate in multi trials (licks per second per trial)')
xlim([-2.5 1])
ylim([0 6])

suptitle(strcat(sessionData.Mouse,'-',sessionData.Date))

%% How to compare psychometric curves ? Place the 50 saliency point in the middle (translation), divide by 0-50 (scaling)


