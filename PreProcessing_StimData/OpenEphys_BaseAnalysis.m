% %% Notes for Claude (and you, dear reader)
% 
% INPUTS:
%   data - OpenEphys recordings in .Rhythm Data format. (a recording is
%          synonomous with a consecutive series of trials for a given condition.
%          Sometimes you will start/stop recording without changing the
%          file name or condition and this creates subfolders in that recording)
%   
%   animal - PROVIDED BY USER. string providing animal name (recommended naming format
%            is YYYYMMDD-p# [dateOfRecording-postnatalDay]
%   
%   stim - PROVIDED BY USER. string providing name of the stimulus in this recording. 
% 
%   other assorted parameters set by the user
% 
% OUTPUTS:
%    [Note: ProbeInfo has the prefix animal-]
%     
%    ProbeInfo.mat - struct containing several fields with information
%                    relevant to the experiment, recording details, and plotting of data.
%                    It is initialized in this script and can be further updated later. Has
%                    the following fields:
%                      
%                    .Animal = animal; %animal name used for file labels and figures
%                    .Areas = areas; %which brain areas are associated with each probe
%                    .ProbeNum = probenum; %number of probes in the recording
%                    .ProbeMaps = probemaps; %2d matrices reflecting shape of each probe and associated channel IDs in data (IDs correspond to rows of raw data)
%                    .ChanIds = chan_ids; %vector of channel IDs in order
%                    .Chans = chans; %number of recording channels
%                    .TTLch = TTLch; %channel ID of TTL trigger channel. should be chans+1 unless something wierd in the GUI during recording
%                    .poi = poi; %which probes in the recording are to be analyzed?
%                    .ProbeIds = probeids; %channel IDs within each probe (starts at 1 : number of channels on probe)
%                    .Ch_Remove = {}; %empty field to be modified later using OpenEphys_editProbeInfo_ChRemove in case some channels are not useful (e.g. outside brain/broken) 
%  
%     [Note: all following output files have the prefix animal-stim-]
%     LFP.mat - mat file containing: 
%               stim_lfp_stimchunks - 3D array of lfp data split into
%                                     trials. Data is downsampled to 1kHz
%                                     from original 30kHz in the recording.
%                                     size[trials x trial length x channels] 
%                  stim_times - vector of stimulus onset timepoints from raw
%                               data.  size[trials]
%                  tr_keep  -  vector of trials considered "good" according to 
%                              user. By default, this contains all trials and
%                              is modified later from a different script
%                              size[trials]
%                  tr_remove - the complement to tr_keep. By default this is
%                              left empty and modified later by a differen
%                              script size[empty]
% 
%    CSD_results.mat - mat file containing:
%                        stim_CSD - cell array containing cell arrays with
%                                   mean CSD across trials (derived from LFP data)
%                                   for each shank in each
%                                   probe. size{1 x probe number}{number of
%                                   shanks}[trial length x channels per shank]
%                        tr_keep  - vector of trials considered "good" according to 
%                                   user. By default, this contains all trials and
%                                   is modified later from a different script
%                                   size[trials]
%                        tr_remove - the complement to tr_keep. By default this is
%                                    left empty and modified later by a differen
%                                    script size[empty]
%                        
% 
%    TF_results.mat - mat file containing:
%                       stim_tf - 4D array containing complex-valued output
%                                 of morlet wavelet convolution on LFP data
%                                 data. size[channels x frequencies x time
%                                 x trials]
%                       tr_keep  - vector of trials considered "good" according to 
%                                  user. By default, this contains all trials and
%                                  is modified later from a different script
%                                  size[trials]
%                       tr_remove - the complement to tr_keep. By default this is
%                                   left empty and modified later by a differen
%                                   script size[empty]
% 
%   spiking_results.mat - mat file containing:
%                           stim_spike_stimchunks - 3D array of MUA data,
%                                                   where each timepoint is labelled as 0 or 1, for
%                                                   no spike or spike. Data is in original 30kHz
%                                                   sampling rate.
%                                                   size[trials x 30kHz trial length x channels]
%                                                 
%                           tr_keep  - vector of trials considered "good" according to 
%                                      user. By default, this contains all trials and
%                                      is modified later from a different script
%                                      size[trials]
%                           tr_remove - the complement to tr_keep. By default this is
%                                       left empty and modified later by a differen
%                                       script size[empty]
% 
%    
%    [Note: all output figures also have brain area added to prefix
%    animal-stim-area(i)]
%    
%    Probe-lfp_results - figure of mean LFP for each channel, plotted in
%                        shape of probe (ex. 8x8 probe has an 8x8
%                        arrangement of subplots). There is one per probe
%    Probe-CSD_results - figure of mean CSD for each shank, plotted in
%                        shape of probe (ex. 8x8 probe has an 1x8
%                        arrangement of subplots). There is one per probe
%    Probe-TF_results -  figure of mean frequency power in time domain 
%                        for each channel, plotted in shape of probe 
%                        (ex. 8x8 probe has an 8x8 arrangement of
%                        subplots). There is one per probe
%    Probe-Spikemean_results -  figure of mean MUA spike rate for each channel, 
%                               plotted in shape of probe (ex. 8x8 probe has an 8x8 arrangement of
%                               subplots). There is one per probe
%    Probe-SpikeRaster_results -  figure of MUA raster for all trials for each channel, 
%                                 plotted in shape of probe (ex. 8x8 probe has an 8x8 arrangement of
%                                 subplots). There is one per probe
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 

%% NEEDS INPUT FIRST
clear all
close all
tic
%% NEEDS INPUT FIRST:  SET UP ANALYSIS INFO
%where is the data stored?
load_directory = 'E:\Roy\Silicon Probe Raw Data\20260616-p7'; %recommend to choose folder for whole experiment

%where will the output be saved?
save_directory = 'E:\Roy\Processed Silicon Probe Data'; %this should stay the same for all animals/recordings

%clear saved data/variables? (likely necessary for memory)
mr_clean = 1;


%%%Analysis details%%%

%which analyses to run?
LFP = 1; %local field potential
CSD = 1; %current source density
TF = 1; %time frequency
MUA = 1; %Spiking

%want to check for bad trials?
TrRemove = 0;

%Want to check for bad channels/channels outside brain? (will not change
%plotting though)
%NOTE: IF YOU ALREADY HAVE PROBEINFO FOR THIS ANIMAL, YOU MAY ALREADY HAVE
%THIS DONE AND CAN KEEP IT 0;
ChRemove = 0;

%Want to manually adjust TimeFrequency Analysis settings?
TFmyfault = 0; 
    % 0 = default settings
    % 1 = manually adjust settings
    % 
    %Allows adjustments to frequency range, wavelet cycle numbers, and
    %baseline period

%want to manually adjust MUA Analysis settings?
MUAmyfault = 1; 
    % 0 = default settings
    % 1 = manually adjust settings
    % 
    %Allows adjustments to threshold setting


%generate summary plots?
showme = 1;

%save summary plots automatically? (showme must be 1 for this to work)
plotsave = 1; 
    %you can still manually save and create your own filename. 
    % Doing the automatic version will create one figure file per animal-condition
    %in other words, it will overwrite previously autosaved figures of the same
    %condition, analysis, and animal


%want to edit plotting settings?
artistemode = 0;
    % 0 = default settings
    % 1 = manually adjust settings
    % 
    %Allows adjustments to plotting window, plot title, axis labels, etc.


%will you keep all trials?
keepem = 1;

%Do you want to have input for certain settings?
%default or myfault? %%% require input here later. for now just ignore it

%% NEEDS INPUT FIRST:  SET UP RECORDING INFO
%%%Experiment details%%%
animal = '20260616-p7'; %animal id (suggested format YYYYMMDD-p#, where p# is age). This will be the beginning of all file names
stim = 'whisker'; %whatever you call this stim. It will be part of all file names so choose wisely


reorder = 0; %need to reorder the channels? (only if there's a mistake in GUI mapping)

%for recordings done with wrong mapping (uncomment the one you need)
    %load 'pre202406_64_32_channel_re_index'
    %load '202408_64_64_channel_re_index_v5.mat'
    %load 20240708_64_32_channel_re_index.mat


%use already existing probe information for this animal?
Cogito_ProbeSum = 0; %this is not automated in case you want to make a new one
    
    %RUN THE LINE BELOW TO DOUBLE CHECK THAT THE PROBE INFO EXISTS!!!!
    %isfile(fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']))



%how many probes? (Fill out regardless of Cogito_ProbeSum)
probenum = 2;

%what probe layout? (Use the key below. NOTE: order matters here) 
probevariant = 6; 
    %single probe experiments
        % 0 = 1x16 %%%not a real option yet
        % 1 = 4x8
        % 2 = 8x8
    %double probe experiments
        % 3 = 4x8 + 4x8
        % 4 = 4x8 + 8x8
        % 5 = 8x8 + 4x8
        % 6 = 8x8 + 8x8
        % 7 = 1x16 + 4x8 %%%not a real option yet
        % 8 = 4x8 + 1x16 %%%not a real option yet
        % 9 = 1x16 + 8x8 %%%not a real option yet
        % 10 = 8x8 + 1x16 %%%not a real option yet


     
%what locations are being recorded? (ONE PER PROBE, IN ORDER)
areas = {'V1','S1'}; %shorthand versions recommended. This will be in file names

%would you like to manually add extra labels to each probe?
%extrabels = 0;  %doesn't exist yet but will add at some point maybe

%how many trigger input channels? 
TTL = 1;



%% 


%%%nothing else to fill in%%%












%% Check for/Create Save Directories

%create LFP folder
if LFP == 1
     if isfolder([save_directory '\LFP\' animal]) == false
        mkdir([save_directory '\LFP\' animal ])
    end
end

%create CSD folder
if CSD == 1
    if isfolder([save_directory '\CSD\' animal]) == false
        mkdir([save_directory '\CSD\' animal])  
    end
end

%create TF folder
if TF == 1
    if isfolder([save_directory '\TF\' animal]) == false
        mkdir([save_directory '\TF\' animal])  
    end
end

%create spiking folder
if MUA == 1
    if isfolder([save_directory '\Spiking\' animal]) == false
        mkdir([save_directory '\Spiking\' animal])  
    end
end

%create probeinfo folder
if isfolder([save_directory '\ProbeInfo']) == false
        mkdir([save_directory '\ProbeInfo'])  
end

%create figures folder for animal
if isfolder([save_directory '\AnimalFigures\' animal]) == false
    mkdir([save_directory '\AnimalFigures\' animal])  
end
%% get directory (user input)
%%%recording info%%%
%where's the data stored? (find file that has folders named "Record node #")
directory = uipickfiles('FilterSpec',load_directory,'Prompt', 'Choose Folder Containing Recording to Analyze [IMMEDIATE SUBFOLDERS SHOULD BE -Record Node ###-]'); %loads file names'E:\Roy\Ungrouped Whiskers\20251011-p11\20251011-p11-light4l22_40'; %for example

%%
if Cogito_ProbeSum == 0
%%
disp('Making Probes')    
chanmap88 = reshape(1:64, 8, []); %for 8x8 !!DONT CHANGE!!
chanmap48 = reshape(1:32, 8, []);  %for 4x8 !!DONT CHANGE!!
switch probevariant
    case 1
        probemaps = {chanmap48};
    case 2 
        probemaps = {chanmap88};
    case 3 
        probemaps = {chanmap48, chanmap48+32};
    case 4
        probemaps = {chanmap48, chanmap88+32};
    case 5 
        probemaps = {chanmap88, chanmap48+64};
    case 6 
        probemaps = {chanmap88, chanmap88+64};
end

%check that length of probemaps matches number of stated probes


for i = 1:size(probemaps,2)
    chtot(i) = numel(probemaps{i});
end
channels = [1:sum(chtot)];
chan_ids = channels;
chans = length(chan_ids);
TTLch = chans+1;
%which of the probes do you want to look at?
poi = [1:probenum]; %indices of probes of interest



%create an id for each channel by probe
for i = 1:size(probemaps,2)
    probeids{i} = reshape(probemaps{i},[],1);
end

%put info together in a struct
ProbeInfo.Animal = animal; %animal name used for file labels and figures
ProbeInfo.Areas = areas; %which brain areas are associated with each probe
ProbeInfo.ProbeNum = probenum; %number of probes in the recording
ProbeInfo.ProbeMaps = probemaps; %2d matrices reflecting shape of each probe and associated channel IDs in data (IDs correspond to rows of raw data)
ProbeInfo.ChanIds = chan_ids; %vector of channel IDs in order
ProbeInfo.Chans = chans; %number of recording channels
ProbeInfo.TTLch = TTLch; %channel ID of TTL trigger channel. should be chans+1 unless something wierd in the GUI during recording
ProbeInfo.poi = poi; %which probes in the recording are to be analyzed?
ProbeInfo.ProbeIds = probeids; %channel IDs within each probe (starts at 1 : number of channels on probe)
ProbeInfo.Ch_Remove = {}; %empty field to be modified later using OpenEphys_editProbeInfo_ChRemove in case some channels are not useful (e.g. outside brain/broken) 

% SAVE PROBE INFO

fname = sprintf([animal '-ProbeInfo']);
save(fullfile([save_directory '\ProbeInfo\' fname]), 'ProbeInfo'); 

% or load existing probe info
else
    load(fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']))
    chans = ProbeInfo.Chans;
end

%% probe stuff for plotting later
%poiidx = cell2mat(ProbeInfo.ProbeIds(ProbeInfo.poi));% >0;
%poiprobeids = reshape(cell2mat(ProbeInfo.ProbeMaps(ProbeInfo.poi)),1,[]);

%some variables needed to plot data in shape matching probe layout
for i = 1:length(ProbeInfo.poi)
    poilayout{i} = ProbeInfo.ProbeMaps{ProbeInfo.poi(i)};
    poilayoutT{i} = reshape(poilayout{i}',1,[]);
    poimapT{i} = cell2mat(ProbeInfo.ProbeMaps(ProbeInfo.poi(i)))';
end


%%





%% Read in LFP data
if LFP == 1 || CSD == 1 || TF == 1 || TrRemove == 1 || MUAmyfault == 1 || ChRemove == 1
disp('Loading LFP Data')
clear data
session = Session(directory{1});

%Get all pieces of recording and concatenate (if recording was paused then
%unpaused, you can have multiple files for a single recording)
for i = 1:size(session.recordNodes{1,1}. recordings,2)
    datums = session.recordNodes{1,1}. recordings{1,i}.continuous('Acquisition_Board-100.acquisition_board').samples;
    %datums = session.recordNodes{1,1}. recordings{1,i}.continuous('Acquisition_Board-100.Rhythm Data').samples;
    datums = downsample(datums',30)'; %to 1kH, each timepoint is 1 ms
    datums = double(datums).*0.1950; %converts int16 to uV
    data{i} = datums;
    
end
data = horzcat(data{:}); %shape is total channels (normal + TTL) x timepoints
clear datums

disp('Done')
%% set channel number for adaptation to different 

%% Find Stim Timepoints: check signal and get stimulus times
disp('Finding Stim Timepoints')
plot(data(ProbeInfo.TTLch,:)'); %check how the signal looks
normbineTTL = round(data(ProbeInfo.TTLch,:)'./max(data(ProbeInfo.TTLch,:)'));
stim_times = find(normbineTTL == 1); %get stimulus times
%plot(data(129,:)'); %check how the signal looks  2x 64 channel recording
%% Find Stim Timepoints: Get indices for onset (this is only here for really weird recordings)
%stimdat = data(97,:); %get the adc data (stimulator trigger-out signal)
%stimdat = data(chans,:); %get the adc data (stimulator trigger-out signal)thresh = 100; %name a threshold value that indicates a stim trigger
%thresh = -1100; %name a threshold value that indicates a stim trigger
%up_thresh = stimdat >thresh; %find all the points above threshold
%up_idx = find(up_thresh); %get indices of those points
%cross_thresh = stimdat(up_idx - 1)< thresh; %now find which ones are preceded by a below threshold value

%stim_times = up_idx(cross_thresh == 1); %now you have stimulus onset times


%% Select peristimulus epochs for analysis

%choose interval size

pre_second = 5; %how many seconds before stim?
post_second = 10; %how many second after stim?

pre = pre_second * 1000; %change if sampling rate changes
post = (post_second * 1000)-1; %ditto
dur = pre+post+1; %trial epoch length
start = pre;

%make an array of time indices to be analyzed

peristim_idx = []; %initialize

for i = 1:length(stim_times) %get them all
    peristim_idx(i,:) = stim_times(i)-pre:stim_times(i)+post;
end

peristim_idx = reshape(peristim_idx,1,[]); %make it into one nice row, epochs end-to-end

%now filter the data columns so we only have the epochs we want

wittled_braindat = data(:, peristim_idx); %matrix of data from selected epochs, with epochs stuck end-to-end (chans x [trials*dur])


%% Split up data into individual trial activity for each channel
disp('LFP: Splitting data into trials')
%create 3D matrix of data

%for i = 1:96
for i = 1:ProbeInfo.Chans
    channel_data = wittled_braindat(i,:);
    geoenginned_ch = reshape(channel_data, size(stim_times,1),dur);
    stim_lfp_stimchunks(:,:,i) = geoenginned_ch; %with dimensions [trials x epoch length x channels]
end

tr_remove = []; %list of trials to remove (recommended not to change at this point)
tr_keep = 1:length(stim_times); %list of trials to keep (should be disjoint from tr_remove). Initialize as all trials
clear geoenginned_ch
%% In case channel order from GUI is wrong, change it here
if reorder == 1
    stim_lfp_stimchunks = stim_lfp_stimchunks(:,:,newchannelidx);
end

end

%% Next, calculate summary data for each channel
if LFP == 1 || CSD == 1 || TrRemove == 1 || MUAmyfault == 1 || ChRemove == 1
%make a variable

stim_lfp_avg = []; %mean of lfp 
stim_lfp_std = []; %standard deviation of lfp


%for Loop, whomever that is.


for i = 1:chans
    stim_lfp_avg(i,:) = mean(stim_lfp_stimchunks(:,:,i),1); % [channels x trial length]
    stim_lfp_std(i,:) = std(stim_lfp_stimchunks(:,:,i),1);  % [channels x trial length]
end 

end


%% 




% The following 3 sections are establishing some extra info based user input
% at the beginning of script and are not parts of the standard LFP
% analysis. 










%% Plot avgs in probe layout to see if needs artifact removal and to select interesting channels
if ChRemove == 1

%start asking questions
ishsokay = 0;
while ishsokay == 0


    %plot the current trial means
    for prb = 1:size(poilayoutT,2)
    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(stim_lfp_avg(poimapT{prb}(chan),start-500:start+1500));
        %xlim();
        title(['CH',num2str(loc)])
        ylim([-200 200]);
        xline(500)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg Response, ALL TRIALS'])
        %title(['CH',num2str(loc)])
        end
        end
    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(stim_lfp_stimchunks(1,start-500:start+1500,poimapT{prb}(chan)));
        %xlim();
        title(['CH',num2str(loc)])
        ylim([-500 500]);
        xline(500)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Single Trial Response'])
        %title(['CH',num2str(loc)])
        end
        end


%get user input for channel choice
dlgtitle = 'You requested Channel Removing option';
promt = {'What channels look bad? (write ch IDs)'};
fieldsize = [1 150];
definput = {''};
opts.Resize = 'on';
opts.WindowStyle = 'normal';
chcheck = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
chBad = str2num(chcheck{1});

clf

%plot user chosen channels mean lfp and trials
 figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        if sum(chBad==loc) > 0
        plot(stim_lfp_avg(poimapT{prb}(chan),start-500:start+1500),'Color','r');
        else
        plot(stim_lfp_avg(poimapT{prb}(chan),start-500:start+1500));
        end
        %xlim();
        title(['CH',num2str(loc)])
        ylim([-200 200]);
        xline(500)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg Response, ALL TRIALS'])
        %title(['CH',num2str(loc)])
        end
        end
    
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        if sum(chBad==loc) > 0
        plot(stim_lfp_stimchunks(1,start-500:start+2500,poimapT{prb}(chan)),'Color','r');
        else
        plot(stim_lfp_stimchunks(1,start-500:start+2500,poimapT{prb}(chan)));
        end
        %xlim();
        title(['CH',num2str(loc)])
        ylim([-200 200]);
        xline(500)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Single Trial Response'])
        %title(['CH',num2str(loc)])
        end
        end

%ask if those channels are okay
quest = 'Happy with removing those channels?';
dlgtitle = 'You requested channel removal option';
btn1 = 'Yes, I like it';
btn2 = 'No, lemme pick again';
answer = questdlg(quest,dlgtitle,btn1,btn2,btn2);
switch answer
    case btn1
        ishsokay = 1;
        ProbeInfo.Ch_Remove{prb} = chBad;
    case btn2
        ishsokay = 0;
end
end

close all
end
end

%data_chunks_poi_avg = squeeze(mean(data_chunks_poi,1));
%%
if MUAmyfault == 1
disp('Asking You Questions. Please answer me.')

%start asking questions
ishsokay = 0;
while ishsokay == 0


    %plot the current trial means
    for prb = 1:size(poilayoutT,2)
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(stim_lfp_avg(poimapT{prb}(chan),start-500:start+2500));
        %xlim();
        title(['CH',num2str(loc)])
        ylim([-200 200]);
        xline(500)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg Response, ALL TRIALS'])
        %title(['CH',num2str(loc)])
        end
        end
    end



%get user input for channel choice
dlgtitle = 'You requested Spike threshold setting option';
promt = {'What CHANNELs NUMBERS(s) to see data of?'};
fieldsize = [1 50];
definput = {''};
opts.Resize = 'on';
opts.WindowStyle = 'normal';
chcheck = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
MUAchcheck = str2num(chcheck{1});

clf

%plot user chosen channels mean lfp 
for i = 1:length(MUAchcheck)
    figure()
    title(['Mean of Ch ' MUAchcheck(i)])
    plot(stim_lfp_avg(MUAchcheck(i),:))
end

%ask if those channels are okay
quest = 'Are these useful channels?';
dlgtitle = 'You requested Spike threshold setting option';
btn1 = 'Yes, I like it';
btn2 = 'No, lemme pick new ones';
answer = questdlg(quest,dlgtitle,btn1,btn2,btn2);
switch answer
    case btn1
        ishsokay = 1;
    case btn2
        ishsokay = 0;
end
end

close all
end


%%
if LFP == 1 || TrRemove == 1
disp('Asking You Questions. Please answer me.')

if TrRemove == 1

   

%start asking questions
startover = 1;
while startover == 1


    %plot the current trial means
    for prb = 1:size(poilayoutT,2)
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(stim_lfp_avg(poimapT{prb}(chan),start-500:start+2500));
        %xlim();
        title(['CH',num2str(loc)])
        ylim([-200 200]);
        xline(500)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg Response, ALL TRIALS'])
        %title(['CH',num2str(loc)])
        end
        end
    end


ishsokay = 0;
while ishsokay == 0

%get user input for channel choice
dlgtitle = 'You requested trial removal';
promt = {'What CHANNEL NUMBERS(s) to see all trials of?'};
fieldsize = [1 50];
definput = {''};
opts.Resize = 'on';
opts.WindowStyle = 'normal'
chcheck = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
chcheck = str2num(chcheck{1});



figure()
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle('Here dem trials you asked fer')
tottrs = size(stim_lfp_stimchunks,1);
for tr=1:size(stim_lfp_stimchunks,1)
    subplot(ceil(tottrs/10),10,tr)
    hold on
    for chi = 1:length(chcheck)
    plot(smooth(stim_lfp_stimchunks(tr,start-500:start+2500,chcheck(chi)),10))
    end
    title('Trial', tr)
    ylim([-500,500])
    xline(500)
    hold off
end

%ask if those channels are okay
quest = 'Are these useful channels?';
dlgtitle = 'You requested trial removal';
btn1 = 'Yes, now let me remove trials';
btn2 = 'No, lemme pick new ones';
answer = questdlg(quest,dlgtitle,btn1,btn2,btn2);
switch answer
    case btn1
        ishsokay = 1;
    case btn2
        ishsokay = 0;
end
end

ishsokayrly = 0;
while ishsokayrly == 0
%get user input for channel
dlgtitle = 'You requested trial removal';
promt = {'What TRIALS would you like to remove?'};
fieldsize = [1 50];
definput = {''};
opts.Resize = 'on';
opts.WindowStyle = 'normal'
tr_remove = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
tr_remove = str2num(tr_remove{1});

stim_lfp_stimchunks_temp = stim_lfp_stimchunks;
stim_lfp_stimchunks_temp(tr_remove,:,:) = [];
for i = 1:chans
    stim_lfp_avg_temp(i,:) = mean(stim_lfp_stimchunks_temp(:,:,i),1);
end 

%plot new avg with trials removed
for prb = 1:size(poilayoutT,2)
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
for chan = 1:numel(poilayoutT{prb})
    loc = poilayoutT{prb}(chan);
    if loc > 0
    subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
    plot(stim_lfp_avg_temp(poimapT{prb}(chan),start-500:start+2500));
    %xlim();
    title(['CH',num2str(loc)])
    ylim([-200 200]);
    xline(500)
    sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg Response, TRIALS REMOVED'])
    %title(['CH',num2str(loc)])
    end
end
end



quest = 'Are you satisfied with these changes?';
dlgtitle = 'You requested trial removal';
btn1 = 'Yes, happy as a clam';
btn2 = 'No, lemme change some things';
btn3 = 'Im a flipflopper. I want to keep all trials';
answer = questdlg(quest,dlgtitle,btn1,btn2,btn3,btn2);
switch answer
    case btn1
        ishsokayrly = 1;
        startover = 0;
        tr_keep(tr_remove) = [];
        stim_lfp_avg = [];
        for i = 1:chans
        stim_lfp_avg(i,:) = mean(stim_lfp_stimchunks(tr_keep,:,i),1);
        end 
        clear stim_lfp_stimchunks_temp

    case btn2
        quest = 'What do you want to redo?';
        dlgtitle = 'You requested trial removal';
        btn4 = 'Choose new CHANNELS';
        btn5 = 'Choose new TRIALS';
        answer2 = questdlg(quest,dlgtitle,btn4,btn5,btn5);
        switch answer2
            case btn4
                ishsokayrly = 1;
                startover = 1;
                close all
            case btn5
                ishsokayrly = 0;
        end
    case btn3
        ishsokayrly = 1;
        startover = 0;
        clear stim_lfp_stimchunks_temp

end
close all
end
end
end

end

%% 
% Now returning to the normal LFP analysis







%% plot results
if showme == 1 & LFP == 1
disp('LFP: Making some pictures for you')

if artistemode == 1
    dlgtitle = 'Please choose your plotting settings for the LFP analysis';
    prompt = {'Time Window start/stop (relative to stim time. Limits are [-5000 10000])',...
        'Set Y-axis Limits (uV)',...
        'Use default figure titles? (0=no,1=yes)',...
    'If above is 0, write new figure title (note title prefix will still be ---Animal_AreaOfProbe---)',...
    'Include xline for stimulus onset? (0=no,1=yes)',...
    'Inclue CH# subplot titles? (0=no,1=yes)'};
    fieldsize = [1 50; 1 50; 1 50; 1 150; 1 50;1 50];
    opts.Resize = 'on';
    opts.WindowStyle = 'normal'
    definput = {'-500 2500', '-500 500', '1', 'Your mother was a hamster and your father smelt of elderberries!', '1', '1'};
    mansettings = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
    
    inputxlimits = str2num(mansettings{1});
    inputylimits = str2num(mansettings{2});
    time = [inputxlimits(1):inputxlimits(2)];
    close all
    %plot the current trial means
    for prb = 1:size(poilayoutT,2)
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(time,stim_lfp_avg(poimapT{prb}(chan),start+inputxlimits(1):start+inputxlimits(2)));
        %xlim();
        title(['CH',num2str(poimapT{prb}(chan))])
        ylim(inputylimits);
        xlim(inputxlimits);
        if str2num(mansettings{4}) == 1
        xline(0)
        else
        end
         if str2num(mansettings{3}) == 0
             sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' mansettings{4}])
         else
             sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg ' stim ' Response'])
         end
         if str2num(mansettings{3}) == 1
            title(['CH',num2str(loc)])
         else
         end
        end
        end
        
        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'lfp_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end
    
    end

else
    time = [-500:2500];
    close all
    %plot the current trial means
    for prb = 1:size(poilayoutT,2)
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(time,stim_lfp_avg(poimapT{prb}(chan),start-500:start+2500));
        %xlim();
        title(['CH',num2str(poimapT{prb}(chan))])
        ylim([-500 500]);
        xlim([-500 2500]);
        xline(0)
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe' ' Avg ' stim ' Response'])
        title(['CH',num2str(loc)])
        end
        end
        
        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'lfp_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end
    
    end
end
end
%% Save LFP Results
if LFP == 1
 disp('LFP: Saving LFP results')

    fname = sprintf([animal '-' stim '-' 'LFP']);
    save([save_directory '\LFP\' animal '\'  fname], 'stim_lfp_stimchunks', 'stim_times', 'tr_remove', 'tr_keep');     
    

end
%%






%Now beginning the CSD analysis section







%% CSD CSD CSD CSD CSD CSD CSD CSD %%


if CSD == 1
disp('CSD: Calculating CSD from LFP avgs ')


%% Configure Data for CSD analysis
for prb = 1:size(ProbeInfo.poi,2)
    shanknum = numel(poilayout{prb})/size(poilayout{prb},1);
    chanpershank = size(poilayout{prb},1);
    for shank = 1:shanknum  
        shankchannel_average{prb}{shank} = stim_lfp_avg(poilayout{prb}(:,shank),:); % [channels on shank x trial length]
    end

end
%% CSD 
%close all
%for i = 1:size(chanmap,2)
%    for j = 1:8
%        shankchannel_smooth{i}(j,:) = shankchannel_average{i}(j,:);
%    end
%end

for prb = 1:size(ProbeInfo.poi,2)
for q = 1:size(poilayout{prb},2)
    figure();
    shankq = shankchannel_average{prb}{q}; %get shank data
    shankq = shankq/1000000; %uV to volts for the function
    
    CSDq = CSDallshank(shankq', 1000, 100E-6,'inverse', 177E-6);
    %drawnow;
    stim_CSD{prb}{q} = CSDq; %cell {1 x number of probes}{1 x number of shanks on probe} [trial length x channels on shank]
    close
end
 
end

%% PLOT CSD

if showme == 1
disp('CSD: Making plots of CSD for your viewing and reviewing pleasure.')

if artistemode == 1
% for setting settings
hf = figure('Name','Colormap Options','Units','normalized'); 
ax(1) = subplot(6,1,1)
title(ax(1),'parula')
colormap(ax(1),flipud(parula))
hCB = colorbar('north');
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(1).Position = [0.15 .85 0.74 0.05];
hCB.Position =   [0.15 .85 0.74 0.05];
%hf.Position(4) = 0.1000;

ax(2) = subplot(6,1,2)
title(ax(2),'turbo')
colormap(ax(2),flipud(turbo))
hCB = colorbar('north');
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(2).Position = [0.15 .7 0.74 0.05];
hCB.Position =   [0.15 .7 0.74 0.05];
%hf.Position(4) = 0.1000;


ax(3) = subplot(6,1,3)
title(ax(3),'jet')
colormap(ax(3),flipud(jet))
hCB = colorbar('north');
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(3).Position = [0.15 .54 0.74 0.05];
hCB.Position =   [0.15 .54 0.74 0.05];
%hf.Position(4) = 0.1000;

ax(4) = subplot(6,1,4)
title(ax(4),'guppy')
colormap(ax(4),slanCM('guppy'))
hCB = colorbar('north');
%set(gca,'Visible',false,'Title',true)
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(4).Position = [0.15 .4 0.74 0.05];
hCB.Position =   [0.15 .4 0.74 0.05];
%hf.Position(4) = 0.1000;

ax(5) = subplot(6,1,5)
title(ax(5),'iceburn')
colormap(ax(5),fliplr(slanCM('iceburn')))
hCB = colorbar('north');
%set(gca,'Visible',false,'Title',true)
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(5).Position = [0.15 .25 0.74 0.05];
hCB.Position = [0.15 .25 0.74 0.05];
%hf.Position(4) = 0.1000;


ax(6) = subplot(6,1,6)
title(ax(6),'redshift')
colormap(ax(6),fliplr(slanCM('redshift')))
hCB = colorbar('north');
%set(gca,'Visible',false,'Title',true)
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(6).Position = [0.15 .1 0.74 0.05];
hCB.Position = [0.15 .1 0.74 0.05];
%hf.Position(4) = 0.1000;

  
    dlgtitle = 'Please choose your plotting settings for the CSD analysis';
    prompt = {'Time Window start/stop in ms (relative to stim time. Limits are [-5000 10000])',...
              'Use default figure titles? (0=no,1=yes)',...
              'If above is 0, write new figure title (note title prefix will still be ---Animal_AreaOfProbe---)',...
              'Include xline for stimulus onset? (0=no,1=yes)',...
              'Please choose your color palette',...
              'LFP or CSD overlay?'};

    fieldsize = [1 50; 1 50; 1 150; 1 50;1 50; 1 50];
    opts.Resize = 'on';
    opts.WindowStyle = 'normal'
    definput = {'-500 2500', '1', 'Your mother was a hamster and your father smelt of elderberries!', '1','guppy','LFP'};
    mansettings = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
    
    inputxlimits = str2num(mansettings{1});
    time = [inputxlimits(1):inputxlimits(2)];
    colchoice = mansettings{5};
    overlay = mansettings{6};
    switch colchoice
        case 'parula'
            cmap = colormap(flipud('parula'));
        case 'turbo'
            cmap = colormap(flipud('turbo'));
        case 'jet'
            cmap = colormap(flipud('jet'));
        case 'guppy'
           cmap = colormap(slanCM('guppy'));
        case 'iceburn'
            cmap = colormap(fliplr(slanCM('iceburn')));
        case 'redshift'
            cmap = colormap(fliplr(slanCM('redshift')));
    end

close

time = [-500:2500];
plotrange = start+inputxlimits(1):start+inputxlimits(2);

for prb = 1:size(ProbeInfo.poi,2)

    clear clims
    for i = 1:size(cell2mat(ProbeInfo.ProbeMaps(prb)),2)
        for j=1:8
            limi = max(abs(stim_CSD{prb}{i}(plotrange,j)));
            limbit = limi/10;
            clims{i}(j,:) = [-limi-limbit limi+limbit];
        end   
    end

csdmap = cell2mat(ProbeInfo.ProbeMaps(prb));

figure()
set(gcf, 'Position', get(0, 'Screensize'));
if str2num(mansettings{2}) == 1
sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' CSD'])
else
    sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' mansettings{3}])
end
    for i = 1:size(cell2mat(ProbeInfo.ProbeMaps(prb)),2)
    
        subplot(1,size(poimapT{1},2),i)
        contour(time,1:size(poilayout{prb},1),stim_CSD{prb}{i}(plotrange,:)',200);
        set(gca,'Color', [0.0, 0.0, 0.0],'YDir', 'reverse','CLim',clims{i}(j,:)) %'YTick', [],
        ylim([.5 8.5])
        xlim([time(1),time(end)]);
        colormap(cmap)
        c = colorbar('southoutside');
        c.Label.String = 'CSD (uA/mm^3)';
        title('Shank ',i)
        %alpha(.7)
       
        for j = 1:8
            hold on
            switch overlay
                case 'LFP'
                plot(time,(-stim_lfp_avg(csdmap(j,i),plotrange)/max(abs(stim_lfp_avg(csdmap(j,i),plotrange)))/2)+(j),...
                Color='white', LineWidth=.5);
                case 'CSD'
                plot(time,(-stim_CSD{prb}{i}(plotrange,j)/max(abs(stim_CSD{prb}{i}(plotrange,j)))/2)+(j),...
                Color='white', LineWidth=.5);
            end
            if str2num(mansettings{4}) == 1
            xline(0,Color='g',LineWidth=2)
            end
            xlim([time(1),time(end)]);
            hold off
        end

    end

        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'CSD_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end
    
end



else

time = [-500:2500];
plotrange = start-500:start+2500;

for prb = 1:size(ProbeInfo.poi,2)

    clear clims
    for i = 1:size(cell2mat(ProbeInfo.ProbeMaps(prb)),2)
        for j=1:8
            limi = max(abs(stim_CSD{prb}{i}(plotrange,j)));
            limbit = limi/10;
            clims{i}(j,:) = [-limi-limbit limi+limbit];
        end   
    end

csdmap = cell2mat(ProbeInfo.ProbeMaps(prb));



figure()
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' CSD'])
    for i = 1:size(cell2mat(ProbeInfo.ProbeMaps(prb)),2)
    
    
        subplot(1,size(poimapT{1},2),i)
        contour(time,1:size(poilayout{prb},1),stim_CSD{prb}{i}(plotrange,:)',200);
        set(gca,'Color', [0.0, 0.0, 0.0],'YDir', 'reverse','CLim',clims{i}(j,:)) %'YTick', [],
        ylim([.5 8.5])
        xlim([time(1),time(end)]);
        colormap(slanCM('guppy'))
        c = colorbar('southoutside');
        c.Label.String = 'CSD (uA/mm^3)';
        title('Shank ',i)
        %alpha(.7)
       
        for j = 1:8
            hold on
            plot(time,(-stim_lfp_avg(csdmap(j,i),plotrange)/max(abs(stim_lfp_avg(csdmap(j,i),plotrange)))/2)+(j),...
            Color='white', LineWidth=.5);
            xline(0,Color='g',LineWidth=2)
            xlim([time(1),time(end)]);
            hold off
        end

    end

        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'CSD_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end
    
end
end

end

%% Save CSD stuff
disp('CSD: Saving CSD results')

fname = sprintf([animal '-' stim '-' 'CSD_results']);
save([save_directory '\CSD\' animal '\' fname], 'stim_CSD', 'tr_remove', 'tr_keep');

clear shankchannel shankchannel_average

end

%% 




%Now starting Time-Frequency Analysis section






%%  Time-Frequency stuff
if TF == 1
disp('TF: Beginning Time-Frequency Analysis')

%% choose data to analyze


data = single(stim_lfp_stimchunks); %reduce to single to save space and time

trialnum = size(data,1);
%trialnum = 30; %do it manually to limit memory use


%% trim data to right trial length
data = data(:,2000:8000,:);


%% setup wavelet parameters
if TFmyfault == 1

    dlgtitle = 'Please choose your settings for the TF analysis';
    prompt = {'Frequency Range ([start# stop#])',...
        'Wavelet Cycle Range (ex. 4-8 gives better time resolution, 8-12 better frequency resolution)',...
    'Window For Setting Baseline (relative to stimulus. limits -3000 to 0)'};
    fieldsize = [1 100; 1 100; 1 100];
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    definput = {'1 150', '6 10', '-2300 -200'};
    mansettings = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);

    freqsetting = str2num(mansettings{1})
    min_freq = freqsetting(1)
    max_freq = freqsetting(2)
    num_frex = (max_freq-min_freq)+2
    frex = linspace(min_freq, max_freq, num_frex);

    range_cycles = str2num(mansettings{2})

    baseline_diff = str2num(mansettings{3})
    times = -3000:1:3000;
    baseline_window = round(size(times,2)/2)+baseline_diff;

else

    %frequency parameters
    min_freq = 1;
    max_freq = 150;
    num_frex = 151;
    frex = linspace(min_freq, max_freq, num_frex);

    %other wavelet parameters
    range_cycles = [6 10];

    baseline_diff = [-2300 -200]; %endpoints for section of trial to use as baseline
    times = -3000:1:3000; %timepoints relative to stim onset to analyze of trials 
    baseline_window = round(size(times,2)/2)+baseline_diff;

end

%% initialize outputs

    
    stim_tf = zeros(chans, num_frex, length(times), trialnum); %complex data
    stim_tfpower = zeros(chans, num_frex, length(times), trialnum); %power
    stim_tfitpc = zeros(chans, num_frex, length(times), trialnum); %phase
    
%% do the FFT
s = logspace(log10(range_cycles(1)), log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
wavtime = -2:1/1000:2;
half_wave = (length(wavtime)-1)/2;

% FFT parameters
nWave = length(wavtime);
nData = 6001 * trialnum; %trial length * trial number, I think
nConv = nWave + nData -1;


% Concatenate trial data


%concatentate trials for all channels
for j = 1:size(data,3)

    chan_dat = reshape(data(:,:,j), trialnum, size(data,2));
    chan_cat = [1];
    for i = 1:trialnum
        chantri_i = chan_dat(i,:);
        chan_cat = cat(2,chan_cat,chantri_i);
        
    end
    data_cat(j,:) = chan_cat(1,2:end);
end

% FFT on all channels of concatenated trial data

%initialize full data variables
data_tf = cell(1,size(data,3));
data_tfphase = cell(1,size(data,3));
data_prefAngle = cell(1,size(data,3));

for j = 1:size(data,3)    %loop through each channel
   
    dataX = fft(data_cat(j,:),nConv);
    
    %initialize the channel variables
    tf = zeros(length(frex),length(times),trialnum,'single');
    itpc = zeros(length(frex),length(times),'single');
    prefAngle = zeros(length(frex),length(times),'single');

    % perform wavelet convolution
    for fi=1:length(frex) %loop over frequencies

        
        wavelet = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2)); %create wavelet and get its FFT
        waveletX = fft(wavelet,nConv);
        waveletX = waveletX ./ max(waveletX);

    
        as = ifft(waveletX .* dataX); %run convolution in one step
        as = as(half_wave+1:end-half_wave);

        
        as = reshape( as, length(times), trialnum ); %reshape back to time X trials

        tf(fi,:,:) = as; %get complex values
        power(fi,:,:) = abs(as).^2; %compute frequency power
        %phase(fi,:,:)   = angle(as); %compute phase
        itpc(fi,:)      = abs( mean( exp(1i*angle(as(:,tr_keep))) ,2)); %compute itpc
        %prefAngle(fi,:) = angle(mean(exp(1i*angle(as)) ,2)); %compute prefAngle
    
    end
   
    %add channel data to full data cells
    stim_tf(j,:,:,:) = tf; %[channels x frequencies x length(times) x trials]
    stim_tfpower(j,:,:,:) = power; %[channels x frequencies x length(times) x trials]
    %stim_tfphase(j,:,:,:) = phase;
    stim_tfitpc(j,:,:,:) = itpc; %[channels x frequencies x length(times) x trials]
    %stim_prefAngle{j} = prefAngle;
disp(['TF: Running TF analysis on channel ' num2str(j)])
end







%if showme == 1
% convert to decibal and calculate avg (just for plotting purposes)
%clear stim_tfpower_dB stim_tfpower_dB_avg

%for ch = 1:size(stim_tfpower,1)
%    stim_tfpower_dB = [];
%    for i = 1:trialnum
        
%    baseline = mean(squeeze(stim_tfpower(ch,:,baseline_window(1):baseline_window(2),i)),2);
%    stim_tfpower_dB(:,:,i) = downsample(downsample(10*log10(bsxfun(@rdivide, squeeze(stim_tfpower(ch,:,:,i)), squeeze(baseline))),2)',2)';
    
%    end
%    stim_tfpower_dB_avg(ch,:,:) = squeeze(mean(stim_tfpower_dB(:,:,tr_keep),3));
%end

%stim_tfpower_dB_avg = downsample()

% calculate mean of power 
%TF Power
%clear stim_tfpower_dB_avg
%for ch = 1:length(ProbeInfo.ChanIds)
%    stim_tfpower_dB_avg(ch,:,:) = squeeze(mean(stim_tfpower_dB(ch,:,:,tr_keep),4));
%end

end
%% plot TF and ITPC in probe layout
if showme == 1
disp('TF: Plotting TF and ITPC results for you :)')   

% convert to decibal and calculate avg (just for plotting purposes)
clear stim_tfpower_dB stim_tfpower_dB_avg
%
clear stim_tfpower_dB stim_tfpower_dB_avg

for ch = 1:size(stim_tfpower,1)
    disp(['TF: Plotting TF-- Converting power to dB for channel ' num2str(ch)])
    stim_tfpower_dB = [];
    for i = 1:trialnum
        
    baseline = mean(squeeze(stim_tfpower(ch,:,baseline_window(1):baseline_window(2),i)),2);
    stim_tfpower_dB(:,:,i) = downsample(downsample(10*log10(bsxfun(@rdivide, squeeze(stim_tfpower(ch,:,:,i)), squeeze(baseline))),2)',2)';
    
    end
    stim_tfpower_dB_avg(ch,:,:) = squeeze(mean(stim_tfpower_dB(:,:,tr_keep),3));
end

if artistemode == 1
% for setting settings
hf = figure('Name','Colormap Options','Units','normalized'); 
ax(1) = subplot(6,1,1)
title(ax(1),'parula')
colormap(ax(1),'parula')
hCB = colorbar('north');
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(1).Position = [0.15 .85 0.74 0.05];
hCB.Position =   [0.15 .85 0.74 0.05];
%hf.Position(4) = 0.1000;

ax(2) = subplot(6,1,2)
title(ax(2),'turbo')
colormap(ax(2),'turbo')
hCB = colorbar('north');
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(2).Position = [0.15 .7 0.74 0.05];
hCB.Position =   [0.15 .7 0.74 0.05];
%hf.Position(4) = 0.1000;


ax(3) = subplot(6,1,3)
title(ax(3),'jet')
colormap(ax(3),'jet')
hCB = colorbar('north');
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(3).Position = [0.15 .54 0.74 0.05];
hCB.Position =   [0.15 .54 0.74 0.05];
%hf.Position(4) = 0.1000;

ax(4) = subplot(6,1,4)
title(ax(4),'guppy')
colormap(ax(4),fliplr(slanCM('guppy')))
hCB = colorbar('north');
%set(gca,'Visible',false,'Title',true)
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(4).Position = [0.15 .4 0.74 0.05];
hCB.Position =   [0.15 .4 0.74 0.05];
%hf.Position(4) = 0.1000;

ax(5) = subplot(6,1,5)
title(ax(5),'iceburn')
colormap(ax(5),slanCM('iceburn'))
hCB = colorbar('north');
%set(gca,'Visible',false,'Title',true)
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(5).Position = [0.15 .25 0.74 0.05];
hCB.Position = [0.15 .25 0.74 0.05];
%hf.Position(4) = 0.1000;


ax(6) = subplot(6,1,6)
title(ax(6),'redshift')
colormap(ax(6),slanCM('redshift'))
hCB = colorbar('north');
%set(gca,'Visible',false,'Title',true)
set(gca, 'XColor', 'none', 'YColor', 'none','Color','none')
ax(6).Position = [0.15 .1 0.74 0.05];
hCB.Position = [0.15 .1 0.74 0.05];
%hf.Position(4) = 0.1000;

  
    dlgtitle = 'Please choose your plotting settings for the TimeFrequency plot';
    prompt = {'Time Window start/stop in ms (relative to stim time. Limits are [-2999 3000])',...
              'Frequency Range to Plot',...
              'Fixed Color Limits for all subplots?',...
              'IF above is 1, please input preferred limits (reflecting range of frequency power)',...
              'Use default figure titles? (0=no,1=yes)',...
              'If above is 0, write new figure title (note title prefix will still be ---Animal_AreaOfProbe---)',...
              'Please choose your color palette',...
              'Include xline for stimulus onset? (0=no,1=yes)',...
              'Inclue CH# subplot titles? (0=no,1=yes)'};
              

    fieldsize = [1 50; 1 50; 1 50; 1 50; 1 50; 1 150; 1 50;1 50; 1 50];
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    definput = {'-2999 3000','1 150', '0','-15 15' '1', 'Your mother was a hamster and your father smelt of elderberries!','iceburn','1','1'};
    mansettings = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
    
    inputxlimits = str2num(mansettings{1});
    time = round(downsample([inputxlimits(1):inputxlimits(2)]/2,2));
    clims = str2num(mansettings{4});
    colchoice = mansettings{7};
    fran = str2num(mansettings{2});
   
    switch colchoice
        case 'parula'
            cmap = colormap('parula');
        case 'turbo'
            cmap = colormap('turbo');
        case 'jet'
            cmap = colormap('jet');
        case 'guppy'
           cmap = colormap(fliplr(slanCM('guppy')));
        case 'iceburn'
            cmap = colormap(slanCM('iceburn'));
        case 'redshift'
            cmap = colormap(slanCM('redshift'));
    end

close


for prb = 1:size(ProbeInfo.poi,2)
figure()
set(gcf, 'Position', get(0, 'Screensize'));
if str2num(mansettings{5})==1
    sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' TF Power (dB)'], 'Interpreter', 'none')
else
    sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' mansettings{6}], 'Interpreter', 'none')
end


for chan = 1:numel(poilayoutT{prb})
    %loc = poilayoutT{prb}(chan);
    %if loc > 0
    subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
    contour(time, downsample(frex,2), squeeze(stim_tfpower_dB_avg(poimapT{prb}(chan),:,1500+time(2):1500+time(end))),50, 'Fill', 'on')
        set(gca,'ydir','normal', 'ylim', fran)% 'XLim',inputxlimits)
        if str2num(mansettings{4}) == 1
        set(gca,'clim',clim)
        end
        colormap(cmap)
        if str2num(mansettings{9}) == 1
        title(['CH',num2str(poimapT{prb}(chan))])
        end
        xlabel('Time (ms)')
        xticks(0:500:time(end))
        if str2num(mansettings{8}) == 1
        xline(0,'LineWidth',1,'Color','r')
        end
    %end
end
 if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'TF_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
 end
end


else

fran = [1 50];
xlim = [-250:1500];
time = [-250:1500];

for prb = 1:size(ProbeInfo.poi,2)
figure()
set(gcf, 'Position', get(0, 'Screensize'));
sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' TF Power (dB)'], 'Interpreter', 'none')


for chan = 1:numel(poilayoutT{prb})
    %loc = poilayoutT{prb}(chan);
    %if loc > 0
    subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
    contour(time, downsample(frex,2), squeeze(stim_tfpower_dB_avg(poimapT{prb}(chan),:,1250:3000)),50, 'Fill', 'on')
        set(gca,'ydir','normal', 'ylim', fran)% 'XLim',xlim)%, 'clim',[-10 10])
        colormap(slanCM('iceburn'))
        %colorbar('southoutside');
        title(['CH',num2str(poimapT{prb}(chan))])
        xlabel('Time (ms)')
        yticks([1:10:50])
        yticklabels([1 10 20 30 40 50])
        xticks([0 500 1000])
        xticklabels([0 1000 2000])
        xline(0,'LineWidth',1,'Color','r')
    %end
end
 if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'TF_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
 end
end
end

disp('Sorry, for space reasons I will be closing all plots now. Please complain to Roy if you really hate this and think hes dumb for doing it.')
disp('(by the way, if you selected the save plots option, then you can just go open them again)')
close all

%%

disp('TF: Saving TF results')
fname = sprintf([animal '-' stim  '_TF_results','.mat']);
save([save_directory '\TF\' animal '\'  fname], 'stim_tf', 'tr_remove','tr_keep', '-v7.3');
clear xlim ylim

end


%% clean up variables for memory

if mr_clean == 1
disp('Deleting a bunch of stuff to avoid memory loss. Nothing scarier than memory loss. ')
disp('The world becomes unfamiliar. You lose yourself. You lose the ones close to you. ')
disp('And they lose you too. And through your days youll be mostly oblivious, ')
disp('unaware of the hardship youre enduring, of the hardship youre causing. ')
disp('Though perhaps you will feel the slow degradation of the enjoyment you used to get from the world ') 
disp('as your memory takes with it the richness and depth of life. ')
disp('And every now and then, itll sneak in. A little moment that hits you from just the right angle ')
disp('and makes you realize something is off, something is lost.')
disp('But then its gone again. ')
disp('And so are you. ')
clear stim_lfp_stimchunks stim_lfp_avg channel_data stim_lfp_std
clear stim_CSD CSDq shankq
clear stim_tfpower stim_tf stim_tfpower_dB stim_tfpower_dB_avg stim_tfphase stim_tfitpc
clear data_cat data_prefAngle data_tf data_tfphase dataX chan_cat chan_dat itpc phase prefAngle tf wavelet waveletX wavtime wittled_braindat
clear as power 
disp('Okay all done with that :)')
end

%%





%MUA analysis begins here











%% SPIKING ANALYSIS













%% read in data
if MUA == 1
disp('MUA: Beginning MUA analysis')
disp('MUA: Loading data')
clear data
for i = 1:size(session.recordNodes{1,1}. recordings,2)
    datums = session.recordNodes{1,1}. recordings{1,i}.continuous('Acquisition_Board-100.acquisition_board').samples;
    %datums = session.recordNodes{1,1}. recordings{1,i}.continuous('Acquisition_Board-100.Rhythm Data').samples;
    %datums = downsample(datums',3)'; %to 10kH, each timepoint is .1 ms
    fs = 30000; 
    datums = double(datums).*0.1950; %converts int16 to uV
    data{i} =datums;
    %data = session.recordNodes{1,2}. recordings{1,1}.continuous('File_Reader-100.Data-A').samples;
    %data = data(1:97,1:20000000); % 7 zeros
end
data = horzcat(data{:}); % [channels x timepoints in recording]
clear datums
%% Find Stim Timepoints: Get indices for onset
%stimdat = data(97,:); %get the adc data (stimulator trigger-out signal)
disp('MUA: Finding stim timepoints')

figure(); plot(data(ProbeInfo.TTLch,:)'); 
normbineTTL = round(data(ProbeInfo.TTLch,:)'./max(data(ProbeInfo.TTLch,:)'));
stim_times = find(diff(normbineTTL)>0);
clear normbineTTL
%% filter and get MAD 
disp('MUA: Filtering out low frequencies')

data_band = bandpass(data(ProbeInfo.ChanIds,:)',[300 6000],fs, 'ImpulseResponse','iir')'; % [channels x timepoints in recording]
clear data




%% Select peristimulus epochs for analysis

%choose interval size

pre_second = 5; %how many seconds before stim?
post_second = 10; %how many second after stim?

pre = pre_second * fs; %change if sampling rate changes
post = (post_second * fs)-1; %ditto


%make an array of time indices to be analyzed

peristim_idx = []; %initialize

for i = 1:length(stim_times) %get them all
    peristim_idx(i,:) = stim_times(i)-pre:stim_times(i)+post;
end

peristim_idx = reshape(peristim_idx,1,[]); %make it into one nice row, epochs end-to-end

%now filter the data columns so we only have the epochs we want

wittled_braindat = data_band(:, peristim_idx); %matrix of data from selected epochs, with epochs stuck end-to-end



%% calc MAD, set threshold, and filter to remove some duplicates

%calculate median absolute deviation of each channel
MAD = -median(abs(data_band-median(data_band,2))/0.6745,2); 

if MUAmyfault == 1
    disp('MUA: Asking you a question. Please answer me.')

    colmap = cool(7);
    for i = 1:size(MUAchcheck,2)
        figure()
        title(['Example MAD Thresholds on Ch' MUAchcheck(i)])
        hold on
        plot(wittled_braindat(MUAchcheck(i),:),'Color','k')
        yline(MAD(MUAchcheck(i))*4,'color',colmap(1,:))
        yline(MAD(MUAchcheck(i))*6,'color',colmap(2,:))
        yline(MAD(MUAchcheck(i))*8,'color',colmap(3,:))
        yline(MAD(MUAchcheck(i))*10,'color',colmap(4,:))
        yline(MAD(MUAchcheck(i))*12,'color',colmap(5,:))
        yline(MAD(MUAchcheck(i))*14,'color',colmap(6,:))
        yline(MAD(MUAchcheck(i))*16,'color',colmap(6,:))
        legend('','MAD*4', 'MAD*6', 'MAD*8', 'MAD*10', 'MAD*12','MAD*14','MAD*16')
        hold off   
    end
    
    dlgtitle = 'You requested to manually set Spike Threshold';
    prompt = {'Please input a MAD multiplier to set the Spike threshold'};
    fieldsize = [1 50];
    definput = {'10'};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    MUAthresh = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
    MUAthresh = str2num(MUAthresh{1});
    disp('Thanks.')
else
    MUAthresh = 10;
end

%get all timepoints over MAD threshold
spikes = wittled_braindat <= MAD*MUAthresh ;

%binarize by treating timepoint of threshold crossing as the spike
%NOTE: this part could be improved to reflect the more standard spike
%sorting approach of spike detection
spikes = cat(2,zeros(chans,1),diff(spikes')'>0); %size of data_band


%% Split up data into individual trial activity for each channel

%create 3D matrix of data
disp('MUA: Splitting data into trials and making summary stuff')

for i = 1:size(wittled_braindat,1)
    channel_data = spikes(i,:);
    banded_data = wittled_braindat(i,:);
    geoenginned_ch = reshape(channel_data, size(stim_times,1),[]);
    geoband_ch = reshape(banded_data, size(stim_times,1),[]);
    stim_spike_stimchunks(:,:,i) = geoenginned_ch; %with dimensions [trials x epoch length x channels]
    %bandpassed_stimchunks(:,:,i) = geoband_ch;
end


%% In case channel order from GUI is wrong, change it here
if reorder == 1 
    stim_spike_stimchunks = stim_spike_stimchunks(:,:,newchannelidx);
end

%% Next, calculate summary data for each channel

%bin spikes per ms
spikes_ms = [];
for ch = 1:size(stim_spike_stimchunks,3)
    
    
    chspikes = squeeze(stim_spike_stimchunks(tr_keep,:,ch));
    
        for batch = 1:size(stim_spike_stimchunks,2)/(fs/1000)
            
            spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
            spikeper(:,batch) = spikebatchi;
        end
        spikes_ms(:,:,ch) = spikeper;
        stim_spike_avgrate(ch,:) = mean(squeeze(movmean(spikeper,50,2)*1000),1);
        
   
end 



%% Plot all channel averages by shank, just to see


%% make raster stuff
clear trials spikeTimes
for ch = 1:size(stim_spike_stimchunks,3)
    isaspike = find(stim_spike_stimchunks(tr_keep,:,ch)');
    trialnums = repmat(1:size(tr_keep,2),size(stim_spike_stimchunks,2),1);
    trials{ch} = trialnums(isaspike)';
    
    spikeTimes{ch} = (seconds(isaspike)/fs);
end
trialStarts = (repmat(1:size(tr_keep,2),1)*size(stim_spike_stimchunks,2))-(post+1);

%make a plot of it all


%% plot raster
if showme == 1
disp('MUA: Making pictures for you')

    if artistemode == 1
    dlgtitle = 'Please choose your plotting settings for the Spike raster analysis';
    prompt = {'Time Window start/stop in ms (relative to stim time. Limits are [-5000 10000])',...
              'Use default figure titles? (0=no,1=yes)',...
              'If above is 0, write new figure title (note title prefix will still be ---Animal_AreaOfProbe---)',...
              'Include xline for stimulus onset? (0=no,1=yes)',...
              'Inclue CH# subplot titles? (0=no,1=yes)'};
    fieldsize = [1 50; 1 50; 1 150; 1 50;1 50];
    opts.Resize = 'on';
    opts.WindowStyle = 'normal'
    definput = {'-500 2500', '1', 'Your mother was a hamster and your father smelt of elderberries!', '1', '1'};
    mansettings = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
    
    inputxlimits = str2num(mansettings{1});
    time = [inputxlimits(1):inputxlimits(2)];


    for prb = 1:size(poilayoutT,2)
        figure()
        if mansettings{2} == 0
            sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' mansettings{3}], 'Interpreter', 'none')
            else
            sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' ' 'Spike by trial' ], 'Interpreter', 'none')
        end
        set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
            loc = poilayoutT{prb}(chan);
            if loc > 0
            subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
            s = spikeRasterPlot(spikeTimes{poimapT{prb}(chan)}, trials{poimapT{prb}(chan)});
            s.AlignmentTimes = (seconds(trialStarts)/fs);
            if mansettings{5}==1
                s.TitleText = ['CH',num2str(poimapT{prb}(chan))];
            end
            ylabel('Trial')
            s.XLimits=seconds([inputxlimits(1)/1000 inputxlimits(2)/1000]);
            end
        end
    
        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'SpikeRaster_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end
    end

    else

    for prb = 1:size(poilayoutT,2)
        figure()
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' ' 'Spikes by trial' ], 'Interpreter', 'none')
        set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
            loc = poilayoutT{prb}(chan);
            if loc > 0
                subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
                s = spikeRasterPlot(spikeTimes{poimapT{prb}(chan)}, trials{poimapT{prb}(chan)});
                s.AlignmentTimes = (seconds(trialStarts)/fs);
                s.TitleText = ['CH',num2str(poimapT{prb}(chan))];
                ylabel('Trial')
                s.XLimits=seconds([-.5 2.5]);
            end
        end
        
        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'SpikeRaster_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end

    end
    end
end
%% plot spike rate
if showme == 1

    if artistemode == 1
        dlgtitle = 'Please choose your plotting settings for the Spike mean analysis';
        prompt = {'Time Window start/stop in ms (relative to stim time. Limits are [-5000 10000])',...
                  'Set Y-axis Limits (spikes/sec)',...
                  'Smoothing Factor (recommend between 1 and 100, but you do you)',...
                  'Use default figure titles? (0=no,1=yes)',...
                  'If above is 0, write new figure title (note title prefix will still be ---Animal_AreaOfProbe---)',...
                  'Include xline for stimulus onset? (0=no,1=yes)',...
                  'Inclue CH# subplot titles? (0=no,1=yes)'};
        fieldsize = [1 50; 1 50; 1 50; 1 50; 1 150; 1 50;1 50];
        opts.Resize = 'on';
        opts.WindowStyle = 'normal'
        definput = {'-500 2500', '0 500', '20' '1', 'Your mother was a hamster and your father smelt of elderberries!', '1', '1'};
        mansettings = inputdlg(prompt,dlgtitle,fieldsize,definput,opts);
        
        inputxlimits = str2num(mansettings{1});
        inputylimits = str2num(mansettings{2});
        time = [inputxlimits(1):inputxlimits(2)];

 
    for prb = 1:size(poilayoutT,2)
        figure()
        if mansettings{4}==0
            sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' ], 'Interpreter', 'none')
        else
            sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' Spike Avg'], 'Interpreter', 'none')
        end
        set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{1},1),size(poimapT{1},2),chan)
        plot(time,smooth(stim_spike_avgrate(poimapT{prb}(chan),start+inputxlimits(1):start+inputxlimits(2)),20));
        xlim([time(1) time(end)]);
        title(['CH',num2str(poimapT{prb}(chan))])
        ylim([inputylimits(1) inputylimits(2)]);
        if mansettings{6}==1
        xline(0)
        end
        if mansettings{7}==1
        title(['CH',num2str(loc)])
        end
        end
        end
         if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'Spikemean_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
         end
    end
   else
   for prb = 1:size(poilayoutT,2)
   figure()
    sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' stim ' Spike Avg'])
    set(gcf, 'Position', get(0, 'Screensize'));
        for chan = 1:numel(poilayoutT{prb})
        loc = poilayoutT{prb}(chan);
        if loc > 0
        subplot(size(poimapT{prb},1),size(poimapT{prb},2),chan)
        plot([-500:2500],smooth(stim_spike_avgrate(poimapT{prb}(chan),start-500:start+2500),25));
        xlim([time(1) time(end)]);
        title(['CH',num2str(poimapT{prb}(chan))])
        ylim([0 500]);
        xline(0)
        title(['CH',num2str(loc)])
        end
        end
        if plotsave == 1
            fname = sprintf([animal '-' stim '-' ProbeInfo.Areas{prb} '_Probe' '-' 'Spikemean_results']);
            savefig(gcf,[save_directory '\AnimalFigures\' animal '\' fname],'compact')
        end
   
    end
end
end

%% SPIKING: save variables
disp('MUA: Saving MUA results')

fname = sprintf([animal '-' stim '-'  'spiking_results','.mat']);
save([save_directory '\Spiking\' animal '\' fname], 'stim_spike_stimchunks', 'tr_remove','tr_keep','-v7.3');


%close all
disp('MUA: Done')

end

toc