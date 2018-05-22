%% Session Viewer -  Raster plot of behavior
figure()
%Parameters
trig = trialData.stimStart; %align with stimulus onset
win = [-3.500 1.6]; %seconds around trigger

init = 1;
numt = 400; %length(trig)
sortby = 'stimside'; %Choose from : 'stimside' 'type' 'respside'

switch sortby
    case 'stimside'
        subplot(1,3,1)
        title('stimLEFT trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)
        ylabel('Trial Number','FontSize',16)

        subplot(1,3,2)
        title('stimRIGHT trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)
        
        subplot(1,3,3)
        title('catch trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)
        
    case 'type'
        subplot(1,4,1)
        title('VISUAL trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)
        ylabel('Trial Number','FontSize',16)

        subplot(1,4,2)
        title('TACTILE trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)
        
        subplot(1,4,3)
        title('MULTI trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)    
        
        subplot(1,4,4)
        title('CATCH trials','FontSize',16)
        patch([win(1) win(2) win(2) win(1)],[0 0 numt+1 numt+1],[0.5 0.5 0.5]);
        set (gca,'Ydir','reverse')
        xlim(win)
        ylim([init numt+1])
        xlabel('Time - Aligned to Stim onset (s)','FontSize',16)             
        %
    case 'respside'
        %
end

for k= [2 52 142 146 163 165 282 323 334]%init:numt
    
    switch sortby
        case 'stimside'
                if trialData.leftCorrect(k)
                    subplot(1,3,1)
                elseif trialData.rightCorrect(k)
                    subplot(1,3,2)
                else 
                    subplot(1,3,3)
                end
        case 'type'
                if trialData.trialType(k)=='V'
                    subplot(1,4,1)
                elseif trialData.trialType(k)=='T'
                    subplot(1,4,2)
                elseif trialData.trialType(k)=='M' 
                    subplot(1,4,3)
                elseif trialData.trialType(k)=='C'
                    subplot(1,4,4)
                else
                    disp('there is an unmathcing trial-type : not CMTV ')
                end
        case 'respside'
            %complete
    end

    % Find in-window stimStart and stimEnd
    ThisSS = trialData.stimStart(k)-trig(k);
    if ThisSS>win(2)
        ThisSS=win(2);
    end
    if ThisSS<win(1)
        ThisSS=win(1);
    end
    ThisSE = trialData.stimEnd(k)-trig(k);
    if ThisSE>win(2)
        ThisSE=win(2);
    end
    if ThisSE<win(1)
        ThisSE=win(1);
    end
    
  
    % Leave out-of-stim periods grey
    
    % Between stimStart and stimEnd : pastel color depending on trial type 
    if trialData.trialType(k)==84 %Tactile:blue-pastel
        patch([ThisSS ThisSE ThisSE ThisSS],[k k k+0.9 k+0.9],[0.4 0.7 1],'EdgeColor','none');
    elseif trialData.trialType(k)==86 %Visual:red-orange
        patch([ThisSS ThisSE ThisSE ThisSS],[k k k+0.9 k+0.9],[1 0.7 0.5],'EdgeColor','none');
    elseif trialData.trialType(k)==77 %Multi: violet  %before:green-lulo:[0.4 1 0.4]
        patch([ThisSS ThisSE ThisSE ThisSS],[k k k+0.9 k+0.9],[0.7 0.2 0.7],'EdgeColor','none');
    elseif trialData.trialType(k)==67 %Catch:dark-grey
        patch([ThisSS ThisSE ThisSE ThisSS],[k k k+0.9 k+0.9],[0.4 0.4 0.4],'EdgeColor','none');
    end
        
    % Add mouse licks : green=R red=L
    if ~isempty(trialData.lickTime{k,1})
     for i=1:length(trialData.lickTime{k,1});
         lick=trialData.lickTime{k,1}(i)-trig(k);
         if strcmp(trialData.lickSide{k,1}(i),'R')
         patch([lick lick+0.05 lick+0.05 lick],[k k k+1 k+1],[0.8 0 0],'EdgeColor','none');
         elseif strcmp(trialData.lickSide{k,1}(i),'L')
         patch([lick lick+0.05 lick+0.05 lick],[k k k+1 k+1],[0 0.8 0],'EdgeColor','none');    
         end
     end 
    end
    % Add rewards
   
end

suptitle(strcat(sessionData.Mouse,sessionData.Date))