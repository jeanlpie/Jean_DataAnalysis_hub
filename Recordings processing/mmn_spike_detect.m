%% options
tic
datafolder = 'D:\cheetahData\J009_recs'; %data source
outputfolder = 'C:\Scratch\Jean DatFiles';%folder to save matlabData.mat file to (creates new subfolder per session)

sessions = {'2018-05-15_11-37-22_J009_rec1st'};

% channels = {'CSC1', 'CSC2', 'CSC3', 'CSC4',  ... 
%     'CSC5', 'CSC6', 'CSC7', 'CSC8',    ...
%     'CSC9', 'CSC10', 'CSC11', 'CSC12', ...
%     'CSC13', 'CSC14', 'CSC15', 'CSC16',...
%     'CSC17','CSC18','CSC19','CSC20',  ...
%     'CSC21','CSC22','CSC23','CSC24',   ...
%     'CSC25','CSC26','CSC27','CSC28',  ...
%     'CSC29','CSC30','CSC31','CSC32'}; ... %probe 1

% {'CSC1_0001', 'CSC2_0001', 'CSC3_0001', 'CSC4_0001',  ...
%     'CSC5_0001', 'CSC6_0001', 'CSC7_0001', 'CSC8_0001',    ...
%     'CSC9_0001', 'CSC10_0001', 'CSC11_0001', 'CSC12_0001', ...
%     'CSC13_0001', 'CSC14_0001', 'CSC15_0001', 'CSC16_0001',...
%     'CSC17_0001','CSC18_0001','CSC19_0001','CSC20_0001',  ...
%     'CSC21_0001','CSC22_0001','CSC23_0001','CSC24_0001',   ...
%     'CSC25_0001','CSC26_0001','CSC27_0001','CSC28_0001',  ...
%     'CSC29_0001','CSC30_0001','CSC31_0001','CSC32_0001'}; %probe 1

%     'CSC33_0001', 'CSC34_0001', 'CSC35_0001', 'CSC36_0001', ...
%     'CSC37_0001', 'CSC38_0001', 'CSC39_0001', 'CSC40_0001', ...
%     'CSC41_0001', 'CSC42_0001', 'CSC43_0001', 'CSC44_0001', ...
%     'CSC45_0001', 'CSC46_0001', 'CSC47_0001', 'CSC48_0001', ...
%     'CSC49_0001', 'CSC50_0001', 'CSC51_0001', 'CSC52_0001', ...
%     'CSC53_0001', 'CSC54_0001', 'CSC55_0001', 'CSC56_0001', ...
%     'CSC57_0001', 'CSC58_0001', 'CSC59_0001', 'CSC60_0001', ...
%     'CSC61_0001', 'CSC62_0001', 'CSC63_0001', 'CSC64_0001'}; % probe 2

channels = {'CSC1', 'CSC2', 'CSC3', 'CSC4',  ...
    'CSC5', 'CSC6', 'CSC7', 'CSC8',    ...
    'CSC9', 'CSC10', 'CSC11', 'CSC12', ...
    'CSC13', 'CSC14', 'CSC15', 'CSC16',...
    'CSC17','CSC18','CSC19','CSC20',  ...
    'CSC21','CSC22','CSC23','CSC24',   ...
    'CSC25','CSC26','CSC27','CSC28',  ...
    'CSC29','CSC30','CSC31','CSC32', ... %probe 1
    'CSC33', 'CSC34', 'CSC35', 'CSC36', ...
    'CSC37', 'CSC38', 'CSC39', 'CSC40', ...
    'CSC41', 'CSC42', 'CSC43', 'CSC44', ...
    'CSC45', 'CSC46', 'CSC47', 'CSC48', ...
    'CSC49', 'CSC50', 'CSC51', 'CSC52', ...
    'CSC53', 'CSC54', 'CSC55', 'CSC56', ...
    'CSC57', 'CSC58', 'CSC59', 'CSC60', ...
    'CSC61', 'CSC62', 'CSC63', 'CSC64'}; % probe 2

nChunkSamp = 512*2000; %Divide data into chunks of 'nChunkSamp' samples to avoid memory issues. Should be a multiple of 512 to avoid problems

chooseThresholds = 0;
doReref = 1;
doFiltering = 1;
doNotch = 0;
freqLow = 300; %frequency for high-pass filtering
freqHigh = 6000; %frequency for low-pass filtering
filtOrd = 4; %order of butterworth filter used

nump=64; %Number of samples per waveform
prePoints = 22; %number of samples before spike peak

%% Main

for iSess = 1:length(sessions)
    
    outputf =fullfile(outputfolder, sessions{iSess});
    mkdir(outputf)
    cd(outputf)
    fullfilename = [outputf,'\matlabData.mat'];
    MatObj = matfile(fullfilename,'Writable',true);
    
    for iShank = 1%:size(channels,1)
        
        %Read in timestamps and Fs for entire session
        dataset = fullfile(datafolder,sessions{iSess}, [channels{iShank,1} '.ncs']);
        [timeStamps, Fs] = Nlx2MatCSC(dataset,[1 0 1 0 0],0,1);
        
        Fs=Fs(1);
        
        app=(0:1E6/Fs:511E6/Fs)';
        app=repmat(app,1,length(timeStamps));
        timeStamps=repmat(timeStamps,512,1);
        timeStamps=timeStamps+app;
        timeStamps=reshape(timeStamps,1,size(timeStamps,1)*size(timeStamps,2));
        
        nSamples = size(timeStamps,2);
        nChannels = size(channels,2);
        
        %Divide session data into chunks of 'nChunkSamp' samples to avoid memory issues
        chunks = [1:nChunkSamp:nSamples,nSamples+1];
        
        spikes_ts=cell(1,nChannels);
        spikes_waveforms=cell(1,nChannels);
        
        h = waitbar(0,'Running Spike Detection...');
        for iChunk = 1:length(chunks)-1
            
            
            waitbar(iChunk/(length(chunks)-1),h);
            
            csc_ts = timeStamps(chunks(iChunk):chunks(iChunk+1)-1);
            csc_data = zeros(nChannels,size(csc_ts,2));
            
            for iChan = 1:nChannels
                
                %Read in samples from current chunk
                datafile = fullfile(datafolder,sessions{iSess}, [channels{iShank,iChan} '.ncs']);
                data=Nlx2MatCSC(datafile,[0 0 0 0 1],0,4,[timeStamps(chunks(iChunk)), timeStamps(chunks(iChunk+1)-1)]);
                data=reshape(data,1,size(data,1)*size(data,2));
                
                %specific to timeskip bug
                data(csc_ts>=7.1747e+09 & csc_ts<=7.1767e+09) = [];
                
                if doNotch
                   data = notch_pie(data,Fs,50,0.0031); %Notch filter at 50Hz 
%                    data = notch_pie(data,Fs,100,0.0031); %Notch filter at 100Hz 
                end
                
                if doFiltering
                    [B,A]=butter(filtOrd,[freqLow freqHigh]/(Fs/2));
                    csc_data(iChan,:)=filtfilt(B,A,data);
                else
                    csc_data(iChan,:) = data;
                end
                
            end
            
            if doReref
                
                data_mean_probe1 = repmat(mean(csc_data(1:32,:)),32,1);
                data_mean_probe2 = repmat(mean(csc_data(33:64,:)),32,1);
                csc_data(1:32,:) = csc_data(1:32,:)-data_mean_probe1;
                csc_data(33:64,:) = csc_data(33:64,:)-data_mean_probe2;

%                 data_mean = repmat(mean(csc_data),nChannels,1);
%                 csc_data = csc_data-data_mean;

            end
            
            clear data data_mean_probe1 data_mean_probe2 
            
            
            %% Threshold selection
            
            if iChunk == 1
                if chooseThresholds
                    
                  for i=1:nChannels
                      plot(csc_ts,csc_data(i,:));
                      grid on;
                      stringa=sprintf('Select threshold for channel %d',i);
                      title(stringa,'FontSize',30);
                      stringa=sprintf('Insert threshold for channel %d: ',i);
                      options.Resize='on';
                      options.WindowStyle='normal';
                      thresholds(i)=str2num(cell2mat(inputdlg(stringa,'Threshold selection',1,{'600'},options)));
                      close;                                            
                  end
                    
                    %Save data
                    save(fullfilename,'thresholds','-v7.3')
                else
                    try
                        load(fullfilename,'thresholds')
                    catch
                        for iChan = 1:nChannels  %!!64
                            thresholds(iChan) = 4*std(csc_data(iChan,:)) ;
                            save(fullfilename,'thresholds','-v7.3')
                        end
                    end
                end
            end
            
            %% Spike extraction
            
            
            for i=1:nChannels
                th=thresholds(i);
                if th>0
                    app=find( csc_data(i,:) > (th) );
                    app1=SplitVec(app,'consecutive');
                    app=[];
                    for j=1:length(app1)
                        [maxvalue, maxind]=max(csc_data(i,app1{j}));
                        app(j)=app1{j}(maxind);
                    end
                elseif th<0
                    app=find(csc_data(i,:)<th);
                    app1=SplitVec(app,'consecutive');
                    app=[];
                    for j=1:length(app1)
                        [minvalue, minind]=min(csc_data(i,app1{j}));
                        app(j)=app1{j}(minind);
                    end
                end
                %Check for spikes too close to the beginning or end of
                %chunk
                if ~isempty(app)
                    while app(1)<nump
                        app(1)=[];
                        if isempty(app)
                            break;
                        end
                    end
                    if ~isempty(app)
                        while (length(csc_ts)-app(end))<nump
                            app(end)=[];
                            if isempty(app)
                                break;
                            end
                        end
                    end
                end
                
                %                 for i=1:nChannels
                spikes_ts=csc_ts(app);
                appSpikes=zeros(length(app),nump);
                for j=1:length(app)
                    appSpikes(j,:)=csc_data(i,app(j)-prePoints:app(j)+nump-prePoints-1);
                end
                spikes_waveforms=appSpikes;
                
                if iChunk == 1
                    MatObj.spikes_waveforms(1,i) = {appSpikes};
                    MatObj.spikes_ts(1,i) = {spikes_ts};
                else
                    tmp = cell2mat(MatObj.spikes_waveforms(1,i));
                    MatObj.spikes_waveforms(1,i) = {[tmp; appSpikes]};
                    tmp = cell2mat(MatObj.spikes_ts(1,i));
                    MatObj.spikes_ts(1,i) = {[tmp, spikes_ts]};
                    clear tmp
                end
                %                 end
                %
                
%                                 for i=1:nChannels
%                                     figure(i);
%                                     app=spikes_waveforms{i};
%                                     plot(mean(app));
%                                     grid on;
%                                     stringa=sprintf('Channel %d',i);
%                                     title(stringa);
%                                 end
%                                 pause;
            end
            
        end
    end
end
close(h)
runtime = toc
%load data and save it again
