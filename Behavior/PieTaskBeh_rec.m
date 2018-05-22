%% PieTaskBeh_REC
%Uses condition number
%modified by antonieta and jean
%% Cut trials
cut_end_trials = 1;
numtrials = length(trialData.trialEnd)-cut_end_trials;
tfields = fieldnames(trialData);
for k=1:length(tfields) %clean trials that didnt happen
        if (isnumeric( trialData.(tfields{k}) ) || iscell( trialData.(tfields{k}) )) && length(trialData.(tfields{k}))>numtrials
            trialData.(tfields{k}) = trialData.(tfields{k})(1:numtrials);
        end
end

%% Tweak stim

%Correct if unimodal trials still have a value for the intensity of the non-shown modality.
%Actually pie_task only saves what will be shown. Zero values are therefore nans.
%Convert nans into zeros:
trialData.stimDeflection(isnan(trialData.stimDeflection))=0;
trialData.stimContrast(isnan(trialData.stimContrast))=0;

% %Contrasts in log ?
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

%% Resp
conds = 1:42; %unique(trialData.conditionnumber);
condsnummax = max(conds);
for k = 1:condsnummax
      TotalTrials(k) = nansum(trialData.conditionnumber==conds(k));  
      HitRight(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      HitLeft(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      NoResp(k)  = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k)) / TotalTrials(k);
      Contrasts(k) = mean(trialData.stimContrast(trialData.conditionnumber==conds(k)));
      Deflections(k) = mean(trialData.stimDeflection(trialData.conditionnumber==conds(k)));
      %for fisher
      RR(k) = nansum(trialData.firstRespLFR=='R' & trialData.conditionnumber==conds(k));
      LL(k) = nansum(trialData.firstRespLFR=='L' & trialData.conditionnumber==conds(k));
      NN(k) = nansum(isnan(trialData.firstRespLFR) & trialData.conditionnumber==conds(k));
end

%% Fisher test
%fisher test
comp = [1,16 ; 2,17 ; 3,18 ; 8,19 ; 12,20 ; 13,21; 32,42];
for k = 1:size(comp,1)
xL = [LL(comp(k,1)),TotalTrials(comp(k,1))-LL(comp(k,1)) ; LL(comp(k,2)) , TotalTrials(comp(k,2))- LL(comp(k,2)) ];
xR = [RR(comp(k,1)),TotalTrials(comp(k,1))-RR(comp(k,1)) ; RR(comp(k,2)) , TotalTrials(comp(k,2))- RR(comp(k,2)) ];
xN = [NN(comp(k,1)),TotalTrials(comp(k,1))-NN(comp(k,1)) ; NN(comp(k,2)) , TotalTrials(comp(k,2))- NN(comp(k,2)) ];
[hL(k),pL(k)] = fishertest(xL);
[hR(k),pR(k)] = fishertest(xR);
[hN(k),pN(k)] = fishertest(xN);
end

%% Figures - Unimodal OR 1D
onlyvis= [7 6 5 1 2 3 4];       onlyvis_opto = [16 17 18];
onlytac = [11 10 1 8 9];        onlytac_opto = [16 19];
facilvis = [15 14 1 12 13];     facilvis_opto = [16 20 21];

%VIS ONLY
figure()
hold on
plot(Contrasts(onlyvis), HitRight(onlyvis),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Contrasts(onlyvis), HitLeft(onlyvis),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Contrasts(onlyvis), NoResp(onlyvis),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
legend({'LickRight','LickLeft','NoResponse'})
%opto
if sum(ismember(onlyvis_opto,unique(trialData.conditionnumber)))>0
plot(Contrasts(onlyvis_opto), HitRight(onlyvis_opto),'Color',[0.8 0 1],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Contrasts(onlyvis_opto), HitLeft(onlyvis_opto),'Color',[0 0.8 1],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Contrasts(onlyvis_opto), NoResp(onlyvis_opto),'Color',[.7 .7 1],'LineWidth',3,'Marker','o','MarkerSize',5)
legend({'LickRight','LickLeft','NoResponse','LickRight_opto','LickLeft_opto','NoResponse_opto'})
end
ylim([0 1])
title('Mouse choices - Visual Only')
ylabel('Percent choice')
xlabel('<-LEFT      Contrast      RIGHT->')


%TAC ONLY
figure()
hold on
plot(Deflections(onlytac), HitRight(onlytac),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Deflections(onlytac), HitLeft(onlytac),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Deflections(onlytac), NoResp(onlytac),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
legend({'LickRight','LickLeft','NoResponse'})
%opto
if sum(ismember(onlytac_opto,unique(trialData.conditionnumber)))>0
plot(Deflections(onlytac_opto), HitRight(onlytac_opto),'Color',[0.8 0 1],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Deflections(onlytac_opto), HitLeft(onlytac_opto),'Color',[0 0.8 1],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Deflections(onlytac_opto), NoResp(onlytac_opto),'Color',[.7 .7 .1],'LineWidth',3,'Marker','o','MarkerSize',5)
legend({'LickRight','LickLeft','NoResponse','LickRight_opto','LickLeft_opto','NoResponse_opto'})
end
ylim([0 1])
title('Mouse choices - Tactile Only')
ylabel('Percent choice')
xlabel('<-LEFT      Deflection      RIGHT->')
legend({'LickRight','LickLeft','NoResponse'})

%FACILVIS
figure()
hold on
plot(Contrasts(facilvis), HitRight(facilvis),'Color',[0.8 0 0],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Contrasts(facilvis), HitLeft(facilvis),'Color',[0 0.8 0],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Contrasts(facilvis), NoResp(facilvis),'Color',[.7 .7 .7],'LineWidth',3,'Marker','o','MarkerSize',5)
legend({'LickRight','LickLeft','NoResponse'})
%opto
if sum(ismember(facilvis_opto,unique(trialData.conditionnumber)))>0
plot(Contrasts(facilvis_opto), HitRight(facilvis_opto),'Color',[0.8 0 1],'LineWidth',3,'Marker','o','MarkerSize',5)
plot(Contrasts(facilvis_opto), HitLeft(facilvis_opto),'Color',[0 0.8 1],'LineWidth',3,'Marker','o','MarkerSize',5) 
plot(Contrasts(facilvis_opto), NoResp(facilvis_opto),'Color',[.7 .7 1],'LineWidth',3,'Marker','o','MarkerSize',5)
legend({'LickRight','LickLeft','NoResponse','LickRight_opto','LickLeft_opto','NoResponse_opto'})
end
ylim([0 1])
title('Mouse choices - Visual Facilitated')
ylabel('Percent choice')
xlabel('<-LEFT      Contrast      RIGHT->')


%% Useful : make conversion matrix to place conditions in 2d.

CondMat = [ 36 35 34 11  0  0  0 ;...
            33 15 14 10 22 23  0 ;...
            7  6  5  1  2  3  4 ;...
            0 25 24  8 12 13 29 ;...
            0  0  0  9 30 31 32 ];    

CondMatOpto = [ 0  0  0  0  0  0  0 ;...
                0  0  0 26 27 28  0 ;...
                0  0  0 16 17 18 37 ;...
                0  0  0 19 20 21 39 ;...
                0  0  0 38 40 41 42 ];
        
SelectedConds   = 1:42; %1:condsnum
HitLeft2D       = nan(5,7);
HitRight2D      = nan(5,7);
NoResp2D        = nan(5,7);
Totr            = nan(5,7);
HitLeft2D_opto       = nan(5,7);
HitRight2D_opto      = nan(5,7);
NoResp2D_opto        = nan(5,7);
TotrOpto             = nan(5,7);
for k = SelectedConds
    HitLeft2D(CondMat==k)   = HitLeft(k);
    HitRight2D(CondMat==k)  = HitRight(k);
    NoResp2D(CondMat==k)    = NoResp(k);
    Totr(CondMat==k)        = TotalTrials(k);
    HitLeft2D_opto(CondMatOpto==k)       = HitLeft(k);
    HitRight2D_opto(CondMatOpto==k)      = HitRight(k);
    NoResp2D_opto(CondMatOpto==k)        = NoResp(k);
    TotrOpto(CondMatOpto==k)             = TotalTrials(k);
end

NVis = size(CondMat,2);
NTac = size(CondMat,1);
posx=[];
for tt=1:NVis
posx = [ posx repmat(tt,1,NTac)];
end
posy = [repmat(1:NTac,1,NVis)];

figure('Units','normalized','pos',[0.05 0.4 0.8 0.5])
colormap(summer)
subplot(2,3,1)
imagesc(HitLeft2D*100,'AlphaData',~isnan(HitLeft2D))
xlabel('VIS')
ylabel('TAC')
text(posx(~isnan(HitLeft2D)),posy(~isnan(HitLeft2D)),num2str(int8(round(( HitLeft2D(~isnan(HitLeft2D(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
title('Percent choice Left')
set(gca,'clim',[0,100],'YTick',1:length(Deflections([9 8 1 10 11])),'YTickLabels',Deflections([9 8 1 10 11]),'XTick',1:length(Contrasts([7 6 5 1 2 3 4])),'XTickLabels',Contrasts([7 6 5 1 2 3 4]));


subplot(2,3,2)
imagesc(HitRight2D*100,'AlphaData',~isnan(HitRight2D))
xlabel('VIS')
ylabel('TAC')
text(posx(~isnan(HitRight2D)),posy(~isnan(HitRight2D)),num2str(int8(round(( HitRight2D(~isnan(HitRight2D(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
set(gca,'clim',[0,100],'YTick',1:length(Deflections([9 8 1 10 11])),'YTickLabels',Deflections([9 8 1 10 11]),'XTick',1:length(Contrasts([7 6 5 1 2 3 4])),'XTickLabels',Contrasts([7 6 5 1 2 3 4]));
title('Percent choice Right')

subplot(2,3,3)
imagesc(NoResp2D*100,'AlphaData',~isnan(NoResp2D))
xlabel('VIS')
ylabel('TAC')
text(posx(~isnan(NoResp2D)),posy(~isnan(NoResp2D)),num2str(int8(round(( NoResp2D(~isnan(NoResp2D(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
set(gca,'clim',[0,100],'YTick',1:length(Deflections([9 8 1 10 11])),'YTickLabels',Deflections([9 8 1 10 11]),'XTick',1:length(Contrasts([7 6 5 1 2 3 4])),'XTickLabels',Contrasts([7 6 5 1 2 3 4]));
title('Percent No Choice')

subplot(2,3,4)
imagesc(HitLeft2D_opto*100,'AlphaData',~isnan(HitLeft2D_opto))
xlabel('VIS')
ylabel('TAC')
text(posx(~isnan(HitLeft2D_opto)),posy(~isnan(HitLeft2D_opto)),num2str(int8(round(( HitLeft2D_opto(~isnan(HitLeft2D_opto(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
set(gca,'clim',[0,100],'YTick',1:length(Deflections([9 8 1 10 11])),'YTickLabels',Deflections([9 8 1 10 11]),'XTick',1:length(Contrasts([7 6 5 1 2 3 4])),'XTickLabels',Contrasts([7 6 5 1 2 3 4]));
title('Percent choice Left Opto')

subplot(2,3,5)
imagesc(HitRight2D_opto*100,'AlphaData',~isnan(HitRight2D_opto))
xlabel('VIS')
ylabel('TAC')
text(posx(~isnan(HitRight2D_opto)),posy(~isnan(HitRight2D_opto)),num2str(int8(round(( HitRight2D_opto(~isnan(HitRight2D_opto(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
set(gca,'clim',[0,100],'YTick',1:length(Deflections([9 8 1 10 11])),'YTickLabels',Deflections([9 8 1 10 11]),'XTick',1:length(Contrasts([7 6 5 1 2 3 4])),'XTickLabels',Contrasts([7 6 5 1 2 3 4]));
title('Percent choice Right Opto')

subplot(2,3,6)
imagesc(NoResp2D_opto*100,'AlphaData',~isnan(NoResp2D_opto))
xlabel('VIS')
ylabel('TAC')
text(posx(~isnan(NoResp2D_opto)),posy(~isnan(NoResp2D_opto)),num2str(int8(round(( NoResp2D_opto(~isnan(NoResp2D_opto(:))) ),2)*100)),'Color',[0 0.15 0.3],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
set(gca,'clim',[0,100],'YTick',1:length(Deflections([9 8 1 10 11])),'YTickLabels',Deflections([9 8 1 10 11]),'XTick',1:length(Contrasts([7 6 5 1 2 3 4])),'XTickLabels',Contrasts([7 6 5 1 2 3 4]));
title('Percent No Choice Opto')

%%
figure() 
colormap(hot);
subplot(2,1,1)
imagesc(Totr);set(gca,'clim',[0,max(Totr(:))]); 
posx = [repmat(1,1,5),repmat(2,1,5),repmat(3,1,5),repmat(4,1,5),repmat(5,1,5),repmat(6,1,5),repmat(7,1,5)];
posy = [repmat(1:5,1,7)];
text(posx,posy,num2str(round(Totr(:))),'Color',[0 0 0.5],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
title('Trial numbers')

subplot(2,1,2)
imagesc(TotrOpto);  set(gca,'clim',[0,max(Totr(:))]);
text(posx,posy,num2str(round(TotrOpto(:))),'Color',[0 0 0.5],'FontSize',11,'FontWeight','bold','HorizontalAlignment','center')
title('Trial numbers OPTO')

%% Difference
figure()
colormap(hot)

%Position of significant comparisons :
[Lj,Li]=find(ismember(CondMat,comp(hL,1)));
[Rj,Ri]=find(ismember(CondMat,comp(hR,1)));
[Nj,Ni]=find(ismember(CondMat,comp(hN,1)));

pvalR = split(num2str(pR(hR),1));
pvalL = split(num2str(pL(hL),1));
pvalN = split(num2str(pN(hN),1));

subplot(1,3,1)
imagesc(abs(HitLeft2D_opto - HitLeft2D),'AlphaData',~isnan(abs(HitLeft2D_opto - HitLeft2D)))
set(gca,'clim',[0,1],'YTick',1:length(unique(Deflections)),'YTickLabels',unique(Deflections),'XTick',1:length(unique(fix(Contrasts*100)/100)),'XTickLabels',unique(fix(Contrasts*100)/100));
title('Left - LeftOpto')
text(Li,Lj,pvalL,'FontSize',16,'Color',[0.6 1 0.6],'HorizontalAlignment','center')

subplot(1,3,2)
imagesc(abs(HitRight2D_opto - HitRight2D),'AlphaData',~isnan(abs(HitRight2D_opto - HitRight2D)))
set(gca,'clim',[0,1],'YTick',1:length(unique(Deflections)),'YTickLabels',unique(Deflections),'XTick',1:length(unique(fix(Contrasts*100)/100)),'XTickLabels',unique(fix(Contrasts*100)/100));
title('Right - RightOpto')
text(Ri',Rj',pvalR,'FontSize',16,'Color',[0.6 1 0.6],'HorizontalAlignment','center')

subplot(1,3,3)
imagesc(abs(NoResp2D_opto - NoResp2D),'AlphaData',~isnan(abs(NoResp2D_opto - NoResp2D)))
set(gca,'clim',[0,1],'YTick',1:length(unique(Deflections)),'YTickLabels',unique(Deflections),'XTick',1:length(unique(fix(Contrasts*100)/100)),'XTickLabels',unique(fix(Contrasts*100)/100));
title('NoResp - NoRespOpto')
text(Ni,Nj,pvalN,'FontSize',16,'Color',[0.6 1 0.6],'HorizontalAlignment','center')

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