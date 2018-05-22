%% Get events (Don't run if you already have events)
%CHANGE THIS :
FileName='D:\cheetahData\2018-01-24_15-21-06_novid\Events.nev'; %!! CAREFUL!! 
%%%%%
FieldSelectionFlags = [1 0 1 0 0];
[events_ts, events_ttl] = Nlx2MatEV( FileName, FieldSelectionFlags, 0 , 1);
clear FieldSelectionFlags FileName
diff_ttl = diff([0 events_ttl]);

%% Solve events
%just stim
stims = [3:19 [3:19]+16384 [3:19]+32768];
endtrial = [31 16415 32799] ;
LFRend = [23 16407 32791];
events_sel_ttl = zeros(1,length(events_ttl)); %sparse but same indexing as ts
k=1;
while k<=length(events_ttl) 
   [bool] = ismember(events_ttl(k),stims);
   if bool == 1
        events_sel_ttl(k) = events_ttl(k);
        if ~isempty(  find( ismember(events_ttl(k+1:end) , endtrial) , 1 )  ) || ~isempty(  find( ismember(events_ttl(k+1:end) , LFRend) , 1 )  )
            k = k + min(find(ismember(events_ttl(k+1:end),endtrial),1) , find(ismember(events_ttl(k+1:end),LFRend),1));
        else
            k=length(events_ttl);
        end
   else
       k=k+1;
   end
end
%clean endtrial
k=1;
while k<=length(events_ttl) 
   [bool1] = ismember(events_ttl(k),endtrial) ;
   if k>1
       if ~isempty(max(events_ts( ismember(events_ttl(1:k-1),endtrial) )))
            [bool2] = ((events_ts(k) - max(events_ts( ismember(events_ttl(1:k-1),endtrial) )))>1E5);
       else
            [bool2] = 1;
       end
   else
       bool2 = 1;
   end
   bool = bool1*bool2;
   
   if bool == 1
        events_sel_ttl(k) = events_ttl(k);
        if ~isempty(  find( ismember(events_ttl(k+1:end) , stims) , 1 )  ) || ~isempty(  find( ismember(events_ttl(k+1:end) , LFRend) , 1 )  )
            k = k + min( find(ismember(events_ttl(k+1:end),stims),1) , find(ismember(events_ttl(k+1:end),LFRend),1));
        else
            k=length(events_ttl);
        end
   else
       k=k+1;
   end
end
%clean LFRend
k=1;
while k<=length(events_ttl) 
     [bool1] = ismember(events_ttl(k),LFRend) ;
   if k>1
       if ~isempty(max(events_ts( ismember(events_ttl(1:k-1),LFRend) )))
            [bool2] = ((events_ts(k) - max(events_ts( ismember(events_ttl(1:k-1),LFRend) )))>1E5);
       else
            [bool2] = 1;
       end
   else
       bool2 = 1;
   end
   bool = bool1*bool2;
   if bool == 1
        events_sel_ttl(k) = events_ttl(k);
        if ~isempty(  find( ismember(events_ttl(k+1:end) , stims) , 1 )  ) || ~isempty(  find( ismember(events_ttl(k+1:end) , endtrial) , 1 )  )
             k = k + min( find(ismember(events_ttl(k+1:end),endtrial),1) , find(ismember(events_ttl(k+1:end),stims),1)) ;
        else
            k=length(events_ttl);
        end
   else
       k=k+1;
   end
end
% licks
events_sel_ttl(events_ttl==16384)=16384; %left
events_sel_ttl(events_ttl==32768)=32768; %right

% Replace
events_ttl = events_sel_ttl;
clear events_sel_ttl bool endtrial stims k LFRend

%% Create events_str (Don't run if already have events_str)
%CAREFUL ! Use the proper TTL to event-type conversion.
events_str=cell(1,length(events_ttl));
%
events_str(ismember(events_ttl,[1 1+16384 1+32768]))= {'blank'};
events_str(ismember(events_ttl,[2 2+16384 2+32768]))= {'blank_opto'};
events_str(ismember(events_ttl,[3 3+16384 3+32768]))= {'visual_L'};
events_str(ismember(events_ttl,[4 4+16384 4+32768]))= {'visual_L_opto'};
events_str(ismember(events_ttl,[5 5+16384 5+32768]))= {'tactile_L'};
events_str(ismember(events_ttl,[6 6+16384 6+32768]))= {'tactile_L_opto'};
events_str(ismember(events_ttl,[7 7+16384 7+32768]))= {'multi_L'};
events_str(ismember(events_ttl,[8 8+16384 8+32768]))= {'multi_L_opto'};
events_str(ismember(events_ttl,[9 9+16384 9+32768]))= {'catch'};
events_str(ismember(events_ttl,[10 10+16384 10+32768]))= {'catch_opto'};
events_str(ismember(events_ttl,[11 11+16384 11+32768]))= {'visual_R'};
events_str(ismember(events_ttl,[12 12+16384 12+32768]))= {'visual_R_opto'};
events_str(ismember(events_ttl,[13 13+16384 13+32768]))= {'tactile_R'};
events_str(ismember(events_ttl,[15 15+16384 15+32768]))= {'tactile_R_opto'};
events_str(ismember(events_ttl,[16 16+16384 16+32768]))= {'multi_R'};
events_str(ismember(events_ttl,[17 17+16384 17+32768]))= {'multi_R_opto'};
%
events_str(ismember(events_ttl,[23 23+16384 23+32768]))= {'LFRend'};
events_str(ismember(events_ttl,[25 25+16384 25+32768]))= {'checker'};
events_str(ismember(events_ttl,[31 31+16384 31+32768]))= {'endtrial'};

%% Do the lick stuff properly:
events_ttl(diff_ttl==16384) = 16384;
events_ttl(diff_ttl==32768) = 32768;
%
events_str(events_ttl==32768)= {'lickRight'};
events_str(events_ttl==16384)= {'lickLeft'};

%%
%What to use to align events_str to trialData?
%endtrial sounds nice
%TTL=10ms so let's allow a 10ms dephase + any actual dephase

%choose event
e1 = 'LFRend';

%assuming the first endtrial correspond to the first entry of
%trialData.trialEnd with no dephasage, we can calculate dephasage on the rest.
%trialData is in seconds, converted below into microseconds.
%events_ts is in microseconds.
%zerEV_TS in microseconds. the distance between events_ts
zerEV_TS = events_ts(strncmp(e1,events_str,6)); zerEV_TS = zerEV_TS - zerEV_TS(1);
zerE1 = ( trialData.respwinEnd*1E6 - trialData.respwinEnd(1)*1E6 )'; 
%  A = trialData.stimStart(trialData.trialType=='T'); 
%  zerE1 = ( A - A(1) )'*1E6; 


%since events can't have extra invalid trials, we take events size as the real number of trials, 
%assuming all events have been properly recorded. 
%Well shit, there is a firstrespLFR at trial 661 and length(EV_TS) is 612
%There are un-evented trialEnds :O
%shitshitshit

dephasage = zerEV_TS - zerE1(1:length(zerEV_TS)); 
figure()
plot(zerEV_TS)
hold on

plot(zerE1)
title(strcat('dephasage for .',e1))
%tolerated dephasage due to TTL overlap bug is +-10ms : 2E4 us

%FIRST shit in trial 30 : too negative dephas.
%then real crap starts trial 32 when dephas grows A LOT.
figure()
plot(dephasage)

%Results for 20180124
%last good trial: 620
%problem with trial: 30, 175, 495
%nlx event precedes beh event by : 8E5, 3E5, 3.2E5 microseconds
%meaning 800ms, 300ms, 320 ms.
%EndTrial is a good one to align.
%300


%%RECHECK
%OTHERWISE DONT USE TIMING, USE RANK
%RANK not good if different trial numbers...