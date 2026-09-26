% SummaryAnalysis_MultiMouseMUA_modifiedforallChannels.m
%
% Description: Pools peri-stimulus MUA across mice and conditions using every
%   usable channel on one probe per mouse. The user lists only the mice, the
%   region per mouse, and the condition names, then picks files. The probe is
%   the one whose ProbeInfo.Areas matches the mouse's region; its channels
%   minus ProbeInfo.Ch_Remove are analyzed. Each file's mouse is found from
%   its file name. Spikes are binned to 1 ms and smoothed into a firing rate
%   (50 ms moving mean), then per-trial AUC and peak rate are computed for an
%   early (P1) and late (P2) window after the stimulus. Results go into one
%   long-format table (one row per trial x channel) and are written to a csv
%   for R. No plotting (plot sections were removed; to be rebuilt).
%
% Inputs:
%   User-set config (edited at top of script):
%     mice (cell array of strings, 1 x nMice) - animal IDs, exactly as they
%       start the data file names; order becomes the factor order in R
%     region (cell array of strings, 1 x nMice) - brain area per mouse (e.g. "V1");
%       must match exactly one entry of that mouse's ProbeInfo.Areas (selects the probe)
%     condis (cell array of char, 1 x nCond) - condition names; order becomes the factor order in R
%   Files picked per condition via uipickfiles (superCondis_dir, cell 1 x nCond),
%   from E:\Roy\Processed Silicon Probe Data\Spiking\<animal>\<animal>-<stim>-spiking_results.mat,
%   each containing:
%     stim_spike_stimchunks (0/1, trials x samples x channels) - peri-stimulus
%       spike indicator at 30 kHz. Stim onset is 1 ms bin 5001 (stim_onset):
%       OpenEphys_BaseAnalysis*.m epochs trials as stim_times-pre:stim_times+post
%       with pre = 5 s, so stim_times sits at position pre+1. Channel index is assumed
%       to equal the raw channel ID [inferred: data is built from
%       data(ProbeInfo.ChanIds,:), so this holds when ChanIds = 1:Chans]
%     tr_keep (double, 1 x nKeep) - indices of good trials (only these are analyzed)
%   E:\Roy\Processed Silicon Probe Data\ProbeInfo\<animal>-ProbeInfo.mat, per mouse:
%     ProbeInfo.Areas (cell, 1 x nProbes) - region name per probe
%     ProbeInfo.ProbeMaps (cell, 1 x nProbes) - 2D map of raw channel IDs per probe (0 = no site)
%     ProbeInfo.Ch_Remove (cell, 1 x nProbes) - raw channel IDs to exclude per probe;
%       missing/empty = exclude none. Assumed indexed by probe number
%       [inferred: OpenEphys_BaseAnalysis*.m indexes it by position in poi, which is
%       the same only when poi = 1:ProbeNum (the default)]
%
% Derived (not user-set):
%   condiFileMice (cell, 1 x nCond) - per file, index into mice (from the file name)
%   condiFileNum (cell, 1 x nCond) - per file, Nth file of that mouse in that condition (1..n)
%   mouseProbe (double, 1 x nMice) - probe index used for each mouse
%   chanlist (cell, 1 x nMice) - sorted column vector of raw channel IDs analyzed per mouse
%
% Outputs:
%   tableres (table, nRows x 14) - one row per trial x channel. Columns:
%     P1_AUC, P2_AUC, All_AUC (double) - integral of the rate over P1 (75-350 ms post
%       stimulus), P2 (350-3000 ms) and both (75-3000 ms), in spikes/s*ms
%       (P1wind / P2wind / Allwind = stim_onset + those ms)
%     P1_Peak, P2_Peak (double) - rate at the highest local maximum in the P1 / P2
%       window (lower on both sides, window padded 1 bin each side; a flat top counts
%       once, at its center), in spikes/s. 0 if the rate is 0 throughout the window;
%       NaN if there are spikes but no local max (rate only decaying/rising).
%       Uses islocalmax (MATLAB R2017b+)
%     P1_PeakTime, P2_PeakTime (double) - time of that peak, in ms post stimulus
%       (bin - stim_onset); NaN whenever there is no peak (Peak 0 or NaN)
%     Animal_Name (string), Animal_Num (double) - mouse ID and its index in mice
%     Condition_Name (string), Condition_Num (double) - condition name and its index in condis
%     Condition_FileNum (double) - Nth file of this mouse within this condition (1..n)
%     Region (string) - from region
%     Channel_ID (double) - raw channel ID (unique across probes; map to probe
%       location / ProbeInfo labels downstream)
%   <filename>.csv - tableres written to
%     E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\
%     (file name chosen by the user in a dialog. NOTE: the existing-name check
%     omits '.csv', so an existing csv with the same name IS overwritten)
%
% Dependencies: OpenEphys_BaseAnalysis*.m (writes *-spiking_results.mat and
%   <animal>-ProbeInfo.mat); OpenEphys_editProbeInfo_ChRemove.m (sets Ch_Remove);
%   optionally QuickTrialRemove.m / ConditionalTrialRemove_callable.m
%   (update tr_keep). Requires uipickfiles (File Exchange).

%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%names must match the start of the data file names (animal-stim-spiking_results.mat)
mice = {"20260226-p12"};%, "20260402-p12"};

%which brain area is being recorded for each mouse? (same order as mice)
%used to pick the probe: must match that mouse's ProbeInfo.Areas
region = {"V1"};

condis = {'L_4','LW_4', 'L_8','LW_8','L_12','LW_12','L_15','LW_15'}; %list conditions to be included, named as you'd prefer. Note the order
% as they will be treated as factors later

assert(numel(region) == numel(mice), ...
    'region has %d entries but mice has %d (need one region per mouse)', numel(region), numel(mice))
%%
for con = 1:length(condis)
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\Spiking','Prompt', ['Choose ' condis{con} ' files']);
end
%% match files to mice, and get channels for each mouse from ProbeInfo

%which mouse is each file from? (file names start with animal-)
condiFileMice = cell(1,numel(condis));
condiFileNum = cell(1,numel(condis));
for con = 1:numel(condis)
    for file = 1:numel(superCondis_dir{con})
        [~,fname] = fileparts(superCondis_dir{con}{file});
        m = find(cellfun(@(mouse) startsWith(fname, mouse + "-"), mice));
        assert(numel(m) == 1, 'Condition %d (%s): file "%s" matches %d of the listed mice (need exactly 1)', ...
            con, condis{con}, fname, numel(m))
        condiFileMice{con}(file) = m;
        condiFileNum{con}(file) = sum(condiFileMice{con} == m); %Nth file of this mouse in this condition
    end
end

%channels to use: those on the probe for that mouse's region, minus Ch_Remove
mouseProbe = zeros(1,numel(mice));
chanlist = cell(1,numel(mice));
for m = 1:numel(mice)
    pinfo = load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' char(mice{m}) '-ProbeInfo.mat']),'ProbeInfo');
    prb = find(strcmp(pinfo.ProbeInfo.Areas, region{m}));
    assert(numel(prb) == 1, '%s: region "%s" matches %d probes in ProbeInfo.Areas (need exactly 1)', ...
        mice{m}, region{m}, numel(prb))

    probechans = nonzeros(pinfo.ProbeInfo.ProbeMaps{prb});
    chremove = [];
    if numel(pinfo.ProbeInfo.Ch_Remove) >= prb
        chremove = pinfo.ProbeInfo.Ch_Remove{prb};
    end
    if any(~ismember(chremove, probechans))
        warning('%s: some Ch_Remove{%d} channels are not on probe %d; check Ch_Remove indexing', mice{m}, prb, prb)
    end

    mouseProbe(m) = prb;
    chanlist{m} = setdiff(probechans(:), chremove(:)); %sorted column of raw channel IDs
    disp([char(mice{m}) ': probe ' num2str(prb) ' (' char(region{m}) '), ' num2str(numel(chanlist{m})) ...
        ' channels kept, removed: ' mat2str(intersect(probechans(:), chremove(:))')])
end
clear pinfo
%%
%load the data you will be using (and only that data), one file at a time
superCondis = [];
superCondis_rate = [];
%superCondis_avg = [];
%superCondis_std = [];

triallabel_condiname = []; %metadata for later
triallabel_condinum = [];
triallabel_condifile = [];
triallabel_animalnum = [];
triallabel_animalname = [];
triallabel_region = [];
triallabel_chanID = [];

fs = 30000;
stim_onset = 5001; %1 ms bin of stimulus onset (0 ms post stimulus); epochs are stim_times-pre:stim_times+post
P1wind = stim_onset + (75:350);   %ms post stimulus 75-350
P2wind = stim_onset + (350:3000); %350-3000
Allwind = stim_onset + (75:3000); %75-3000

for con = 1:numel(condis)

    p1auc_con = [];
    p2auc_con = [];
    allauc_con = [];
    p1peak_con = [];
    p2peak_con = [];
    p1peaktime_con = [];
    p2peaktime_con = [];

    for file = 1:numel(superCondis_dir{con})
        disp(['Running Condition ', num2str(con), ' (', condis{con}, '), File ', num2str(file)])
        m = condiFileMice{con}(file);
        chs = chanlist{m};
        dat = matfile(superCondis_dir{con}{file});
        trs = dat.tr_keep;
        ntr = length(trs);
        nch = length(chs);

        %rows are added channel by channel, each with all trials: labels follow the same order
        triallabel_condiname = cat(1,triallabel_condiname,repmat(string(condis{con}),ntr*nch,1));
        triallabel_condinum = cat(1,triallabel_condinum,repmat(con,ntr*nch,1));
        triallabel_condifile = cat(1,triallabel_condifile,repmat(condiFileNum{con}(file),ntr*nch,1));
        triallabel_animalnum = cat(1,triallabel_animalnum,repmat(m,ntr*nch,1));
        triallabel_animalname = cat(1,triallabel_animalname,repmat(string(mice{m}),ntr*nch,1));
        triallabel_region = cat(1,triallabel_region,repmat(string(region{m}),ntr*nch,1));
        triallabel_chanID = cat(1,triallabel_chanID,reshape(repmat(chs',ntr,1),[],1));

        superCondis_filetrials = dat.stim_spike_stimchunks(trs,:,chs);
        nbatch = floor(size(superCondis_filetrials,2)/(fs/1000));

        spikes_rate = zeros(ntr,nbatch,nch);
        %spikes_rateavg = [];
        %spikes_ratestd = [];
        for ch = 1:nch

            chspikes = reshape(superCondis_filetrials(:,:,ch),ntr,[]);
            spikeper = zeros(ntr,nbatch);

            for batch = 1:nbatch

                spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
                spikeper(:,batch) = spikebatchi;

            end

            spikes_rate(:,:,ch) = movmean(spikeper,50,2)*1000;
            %spikes_rateavg(ch,:) = mean(squeeze(spikes_rate(:,:,ch)),1);
            %spikes_ratestd(ch,:) = std(squeeze(spikes_rate(:,:,ch)),1);


            auc = cumtrapz(spikes_rate(:,P1wind,ch),2);
            p1auc_con = cat(1,p1auc_con,auc(:,end));
            auc = cumtrapz(spikes_rate(:,P2wind,ch),2);
            p2auc_con = cat(1,p2auc_con,auc(:,end));
            auc = cumtrapz(spikes_rate(:,Allwind,ch),2);
            allauc_con = cat(1,allauc_con,auc(:,end));

            %Peaks: highest local max in the window (lower on both sides), and its time in ms post stimulus.
            %window is padded 1 bin each side so a decay into / rise out of the window isn't a peak.
            %rounded so movmean plateaus are exactly flat, and a plateau's peak is its center.
            %no local max (only decay/rise) -> Peak and PeakTime NaN; no spikes (rate all 0) -> Peak 0, PeakTime NaN
            seg = round(spikes_rate(:,P1wind(1)-1:P1wind(end)+1,ch),6);
            lm = islocalmax(seg,2,'FlatSelection','center');
            lm(:,[1 end]) = false;
            seg(~lm) = -Inf;
            [peaks,tps] = max(seg,[],2);
            peaktime = tps + P1wind(1) - 2 - stim_onset; %padded index -> bin -> ms post stimulus
            nopeak = isinf(peaks);
            allzero = all(spikes_rate(:,P1wind,ch) == 0,2);
            peaks(nopeak) = NaN;  peaks(allzero) = 0;
            peaktime(nopeak | allzero) = NaN;
            p1peak_con = cat(1,p1peak_con,peaks);
            p1peaktime_con = cat(1,p1peaktime_con,peaktime);

            seg = round(spikes_rate(:,P2wind(1)-1:P2wind(end)+1,ch),6);
            lm = islocalmax(seg,2,'FlatSelection','center');
            lm(:,[1 end]) = false;
            seg(~lm) = -Inf;
            [peaks,tps] = max(seg,[],2);
            peaktime = tps + P2wind(1) - 2 - stim_onset; %padded index -> bin -> ms post stimulus
            nopeak = isinf(peaks);
            allzero = all(spikes_rate(:,P2wind,ch) == 0,2);
            peaks(nopeak) = NaN;  peaks(allzero) = 0;
            peaktime(nopeak | allzero) = NaN;
            p2peak_con = cat(1,p2peak_con,peaks);
            p2peaktime_con = cat(1,p2peaktime_con,peaktime);

        end
        clear superCondis_filetrials
    end
    %superCondis_avg(con,:,:) = spikes_rateavg;
    %superCondis_std(con,:,:) = spikes_ratestd;

    p1auc{con}  = p1auc_con;
    p2auc{con}  = p2auc_con;
    allauc{con} = allauc_con;
    p1peak{con} = p1peak_con;
    p2peak{con} = p2peak_con;
    p1peaktime{con} = p1peaktime_con;
    p2peaktime{con} = p2peaktime_con;
end







%% calculate summary data (peak and AUC)


%calculate auc and peak for p1 and p2 of response

%for con = 1:numel(condis)
%    p1auc_con = [];
%    p2auc_con = [];
%    allauc_con = [];
%    p1peak_con = [];
%    p2peak_con = [];

%    for ch = 1:size(chans,1)

        %AUC
%        auc = cumtrapz(superCondis_rate{con}(:,5150:5350,ch),2);
%        p1auc_con = cat(1,p1auc_con,auc(:,end));
%        auc = cumtrapz(superCondis_rate{con}(:,5350:8000,ch),2);
%        p2auc_con = cat(1,p2auc_con,auc(:,end));
%        auc = cumtrapz(superCondis_rate{con}(:,5020:8000,ch),2);
%        allauc_con = cat(1,allauc_con,auc(:,end));

        %Peaks
%        peaks = max(superCondis_rate{con}(:,5150:5350,ch),[],2);
%        p1peak_con = cat(1,p1peak_con,peaks(:,end));
%        peaks = max(superCondis_rate{con}(:,5350:8000,ch),[],2);
%        p2peak_con = cat(1,p2peak_con,peaks(:,end));
        %peaks = max(superCondis_rate{con}(:,5020:8000),[],2);
        %allpeak{con} = peaks;

%    end

%    p1auc{con}  = p1auc_con;
%    p2auc{con}  = p2auc_con;
%    allauc{con} = allauc_con;
%    p1peak{con} = p1peak_con;
%    p2peak{con} = p2peak_con;
%end
%% make a table that can be easily read out in R for stats and plotting
trialdata_p1auc = [];
trialdata_p2auc = [];
trialdata_allauc = [];
trialdata_p1peak = [];
trialdata_p2peak = [];
trialdata_p1peaktime = [];
trialdata_p2peaktime = [];

for con = 1:numel(condis)
   trialdata_p1auc = vertcat(trialdata_p1auc,p1auc{con});
   trialdata_p2auc = vertcat(trialdata_p2auc,p2auc{con});
   trialdata_allauc = vertcat(trialdata_allauc,allauc{con});
   trialdata_p1peak = vertcat(trialdata_p1peak,p1peak{con});
   trialdata_p2peak = vertcat(trialdata_p2peak,p2peak{con});
   trialdata_p1peaktime = vertcat(trialdata_p1peaktime,p1peaktime{con});
   trialdata_p2peaktime = vertcat(trialdata_p2peaktime,p2peaktime{con});
end


tableres = table(trialdata_p1auc,trialdata_p2auc,trialdata_allauc,trialdata_p1peak,trialdata_p2peak,...
      trialdata_p1peaktime,trialdata_p2peaktime,...
      triallabel_animalname, triallabel_animalnum, triallabel_condiname,triallabel_condinum,triallabel_condifile,...
    triallabel_region,triallabel_chanID,...
    'VariableNames', ["P1_AUC","P2_AUC","All_AUC","P1_Peak","P2_Peak","P1_PeakTime","P2_PeakTime",...
    "Animal_Name","Animal_Num","Condition_Name","Condition_Num","Condition_FileNum","Region","Channel_ID"]);

%% save the table as a csv for R

%give the file a name

dlgtitle = 'What file name to use for the data?';
promt = {'Give the file a name (no spaces please)'};
fieldsize = [1 150];
definput = {''};
opts.Resize = 'on';
opts.WindowStyle = 'normal';
filename = inputdlg(promt,dlgtitle,fieldsize,definput,opts);

%does it exist?
exists = isfile(fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\' cell2mat(filename)]));

while exists == 1
    dlgtitle = 'You stupid, dementia-riddled dumbass. You already made a file called that';
    promt = {'Give the file a name (that doesnt already exist this time)'};
    fieldsize = [1 150];
    definput = {''};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    filename = inputdlg(promt,dlgtitle,fieldsize,definput,opts);

    exists = isfile(fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\' cell2mat(filename)]));

end

writetable(tableres,fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\' cell2mat(filename) '.csv']))
disp('Saved. Go to R, traitor.')


%save(fullfile(['E:\Roy\Processed Silicon Probe Data\BundledAnimalData\20260226-p12_wholeprobe_means']), 'superCondis_avg', 'superCondis_std','superCondis_dir')
