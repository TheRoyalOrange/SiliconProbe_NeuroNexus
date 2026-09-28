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
%     tr_conditional_use (cell array of strings, 1 x nLabels) - names of trial labels made with
%       ConditionalTrialRemove.m to add as Excl_<name> columns; empty = none
%   Files picked per condition via uipickfiles (superCondis_dir, cell 1 x nCond),
%   from E:\Roy\Processed Silicon Probe Data\Spiking\<animal>\<animal>-<stim>-spiking_results.mat,
%   each containing:
%     stim_spike_stimchunks (0/1, trials x samples x channels) - peri-stimulus
%       spike indicator at 30 kHz. Stim onset is 1 ms bin 5001 (stim_onset):
%       OpenEphys_BaseAnalysis*.m epochs trials as stim_times-pre:stim_times+post
%       with pre = 5 s, so stim_times sits at position pre+1. Channel index is assumed
%       to equal the raw channel ID [inferred: data is built from
%       data(ProbeInfo.ChanIds,:), so this holds when ChanIds = 1:Chans]
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1 of stim_spike_stimchunks)
%     tr_remove (0/1, 1 x nTrials) - mask over tr_keep (updated format, as written by
%       ConditionalTrialRemove_callable.m); only tr_keep(~tr_remove) are analyzed.
%       Must be the same length as tr_keep, else the script stops with an error
%     tr_remove_conditional (table, columns Name (string) / TrialIdx (double), optional) -
%       trial labels from ConditionalTrialRemove_callable.m; one Name can have many rows;
%       NaN TrialIdx = evaluated, none flagged; Name absent = file not evaluated for it
%   E:\Roy\Processed Silicon Probe Data\ProbeInfo\<animal>-ProbeInfo.mat, per mouse:
%     ProbeInfo.Areas (cell, 1 x nProbes) - region name per probe
%     ProbeInfo.ProbeMaps (cell, 1 x nProbes) - 2D map of raw channel IDs per probe (0 = no site);
%       column s = shank s, row 1 = top of the shank
%     ProbeInfo.ShankLabels (cell, 1 x nProbes, optional) - per probe a struct (or [] if not labelled);
%       each field holding one number per shank (1 x nShanks, ProbeMaps column order, raw values)
%       is offered as a shank activity score (e.g. SponActivity, LightOnly_mixedIntensity)
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
%   shankscore_use (string, 1 x nScores) - shank activity scores chosen in a list dialog at the
%     start (from those found across the mice); if some mice lack one, a dialog offers
%     abort or continue (those mice get NaN in that column)
%
% Outputs:
%   tableres (table, nRows x (17 + nLabels + nScores)) - one row per trial x channel. Measure columns
%   carry "MUA" so tables from the LFP/MUA/TF scripts can be combined. Columns:
%     P1_MUAAUC, P2_MUAAUC, All_MUAAUC (double) - integral of the rate over P1 (75-350 ms post
%       stimulus), P2 (350-3000 ms) and both (75-3000 ms), in spikes/s*ms
%       (P1wind / P2wind / Allwind = stim_onset + those ms)
%     P1_MUAPeak, P2_MUAPeak (double) - rate at the highest local maximum in the P1 / P2
%       window (lower on both sides, window padded 1 bin each side; a flat top counts
%       once, at its center), in spikes/s. 0 if the rate is 0 throughout the window;
%       NaN if there are spikes but no local max (rate only decaying/rising).
%       Uses islocalmax (MATLAB R2017b+)
%     P1_MUAPeakTime, P2_MUAPeakTime (double) - time of that peak, in ms post stimulus
%       (bin - stim_onset); NaN whenever there is no peak (Peak 0 or NaN)
%     Animal_Name (string), Animal_Num (double) - mouse ID and its index in mice
%     Condition_Name (string), Condition_Num (double) - condition name and its index in condis
%     Condition_FileNum (double) - Nth file of this mouse within this condition (1..n)
%     Trial (double) - original trial number (value from tr_keep); with Animal/Condition/
%       FileNum/Channel_ID, identifies a row across the LFP/MUA/TF tables
%     Region (string) - from region
%     Channel_ID (double) - raw channel ID (unique across probes; map to probe
%       location / ProbeInfo labels downstream)
%     Excl_<name> (double, one per tr_conditional_use name) - 1 = trial flagged under that
%       label, 0 = file evaluated and trial not flagged, NaN = file never evaluated for it
%     (added in their own section after the table is built, so they can be redone alone:)
%     Chan_Shank (double) - shank (ProbeMaps column) the channel is on
%     Chan_Depth (double) - position among that shank's kept (not Ch_Remove) channels,
%       1 = top; NaN if the channel has since been added to Ch_Remove
%     ShankRank_<score> (double, one per shankscore_use) - rank of the channel's shank by that
%       ShankLabels score, 1 = highest (ties -> lower shank number first); NaN if missing
%   <filename>.csv - tableres written to
%     E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR\
%     (file name chosen by the user in a dialog. NOTE: the existing-name check
%     omits '.csv', so an existing csv with the same name IS overwritten)
%
% Dependencies: OpenEphys_BaseAnalysis*.m (writes *-spiking_results.mat and
%   <animal>-ProbeInfo.mat); OpenEphys_editProbeInfo_ChRemove.m (sets Ch_Remove);
%   ConditionalTrialRemove.m / ConditionalTrialRemove_callable.m (tr_remove mask format,
%   tr_remove_conditional labels). Requires uipickfiles (File Exchange).

%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%names must match the start of the data file names (animal-stim-spiking_results.mat)
mice = {"20260226-p12"};%, "20260402-p12"};

%which brain area is being recorded for each mouse? (same order as mice)
%used to pick the probe: must match that mouse's ProbeInfo.Areas
region = {"V1"};

condis = {'L_4','LW_4', 'L_8','LW_8','L_12','LW_12','L_15','LW_15'}; %list conditions to be included, named as you'd prefer. Note the order
% as they will be treated as factors later

%trial labels made with ConditionalTrialRemove.m to add as columns (empty = none)
%each becomes a column Excl_<name>: 1 = flagged, 0 = evaluated & not flagged, NaN = file not evaluated
tr_conditional_use = {}; %e.g. {"whisker_twitch_artifact"}

assert(numel(region) == numel(mice), ...
    'region has %d entries but mice has %d (need one region per mouse)', numel(region), numel(mice))
%% choose which shank activity scores to add as ShankRank_<score> columns
%scores = fields of ProbeInfo.ShankLabels{probe} holding one number per shank (ProbeMaps column),
%for the probe matching that mouse's region
scorenames = cell(1,numel(mice));
for m = 1:numel(mice)
    pinfo = load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' char(mice{m}) '-ProbeInfo.mat']),'ProbeInfo');
    prb = find(strcmp(pinfo.ProbeInfo.Areas, region{m}));
    assert(numel(prb) == 1, '%s: region "%s" matches %d probes in ProbeInfo.Areas (need exactly 1)', ...
        mice{m}, region{m}, numel(prb))

    scorenames{m} = strings(1,0);
    SL = [];
    if isfield(pinfo.ProbeInfo,'ShankLabels') && numel(pinfo.ProbeInfo.ShankLabels) >= prb
        SL = pinfo.ProbeInfo.ShankLabels{prb};
    end
    if isstruct(SL)
        nshank = size(pinfo.ProbeInfo.ProbeMaps{prb},2);
        fn = fieldnames(SL);
        for f = 1:numel(fn)
            if isnumeric(SL.(fn{f})) && numel(SL.(fn{f})) == nshank
                scorenames{m}(end+1) = string(fn{f});
            end
        end
    end
end
clear pinfo SL

allscores = unique([scorenames{:}],'stable');
shankscore_use = strings(1,0);
if isempty(allscores)
    disp('No shank activity scores found in ProbeInfo.ShankLabels for these mice - no ShankRank columns')
else
    [sel,ok] = listdlg('ListString',cellstr(allscores),'SelectionMode','multiple', ...
        'PromptString','Shank activity scores to add as ShankRank_ columns:','ListSize',[350 200]);
    if ok
        shankscore_use = allscores(sel);
    end
    disp(['ShankRank columns: ' char(strjoin(shankscore_use,', '))])
end

%do all mice have the chosen scores?
missingmsg = strings(0,1);
for m = 1:numel(mice)
    miss = setdiff(shankscore_use, scorenames{m}, 'stable');
    if ~isempty(miss)
        missingmsg(end+1,1) = string(mice{m}) + " (" + string(region{m}) + "): " + strjoin(miss,", ");
    end
end
if ~isempty(missingmsg)
    warning('Some mice are missing selected shank scores:\n%s', strjoin(missingmsg, newline))
    answer = questdlg(sprintf('Some mice are missing selected shank scores (those columns will be NaN for them):\n\n%s', ...
        strjoin(missingmsg, newline)), 'Missing shank scores', 'Continue (missing = NaN)', 'Abort', 'Abort');
    if ~strcmp(answer,'Continue (missing = NaN)')
        error('Aborted by user: selected shank scores are missing for some mice')
    end
end
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
triallabel_trialnum = [];
triallabel_trcond = zeros(0,numel(tr_conditional_use)); %one column per tr_conditional_use name

%for the trial label summary printed after the loop
trcond_nflagged = zeros(1,numel(tr_conditional_use));
trcond_noteval = repmat({strings(0,1)},1,numel(tr_conditional_use));

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
        [~,fname] = fileparts(superCondis_dir{con}{file});
        filevars = who(dat);

        %good trials: tr_remove is a 0/1 mask over tr_keep
        tr_keep = dat.tr_keep;
        tr_remove = dat.tr_remove;
        assert(numel(tr_remove) == numel(tr_keep) && all(ismember(tr_remove,[0 1])), ...
            '%s: tr_remove must be a 0/1 mask the same length as tr_keep (updated format)', fname)
        trs = tr_keep(~logical(tr_remove));
        trs = trs(:)';
        ntr = length(trs);
        nch = length(chs);

        %trial labels: 1 = flagged under that name, 0 = evaluated & not flagged, NaN = not evaluated for this file
        %(a NaN TrialIdx row means evaluated with none flagged, and never matches a trial)
        if ismember('tr_remove_conditional',filevars)
            tr_remove_conditional = dat.tr_remove_conditional;
        else
            tr_remove_conditional = table('Size',[0,2],'VariableTypes',{'string','double'},'VariableNames',{'Name','TrialIdx'});
        end
        trflags = nan(ntr,numel(tr_conditional_use));
        for k = 1:numel(tr_conditional_use)
            namerows = strcmp(tr_remove_conditional.Name, tr_conditional_use{k});
            if any(namerows)
                trflags(:,k) = ismember(trs, tr_remove_conditional.TrialIdx(namerows));
                trcond_nflagged(k) = trcond_nflagged(k) + sum(trflags(:,k));
            else
                trcond_noteval{k}(end+1,1) = string(fname);
            end
        end

        %rows are added channel by channel, each with all trials: labels follow the same order
        triallabel_condiname = cat(1,triallabel_condiname,repmat(string(condis{con}),ntr*nch,1));
        triallabel_condinum = cat(1,triallabel_condinum,repmat(con,ntr*nch,1));
        triallabel_condifile = cat(1,triallabel_condifile,repmat(condiFileNum{con}(file),ntr*nch,1));
        triallabel_animalnum = cat(1,triallabel_animalnum,repmat(m,ntr*nch,1));
        triallabel_animalname = cat(1,triallabel_animalname,repmat(string(mice{m}),ntr*nch,1));
        triallabel_region = cat(1,triallabel_region,repmat(string(region{m}),ntr*nch,1));
        triallabel_chanID = cat(1,triallabel_chanID,reshape(repmat(chs',ntr,1),[],1));
        triallabel_trialnum = cat(1,triallabel_trialnum,repmat(trs',nch,1));
        triallabel_trcond = cat(1,triallabel_trcond,repmat(trflags,nch,1));

        %matfile can only read evenly spaced ranges, so read the block spanning the good trials
        %and kept channels once, then pick the ones needed from it in memory
        tr_rng = min(trs):max(trs);
        ch_rng = min(chs):max(chs);
        superCondis_filetrials = dat.stim_spike_stimchunks(tr_rng,:,ch_rng);
        superCondis_filetrials = superCondis_filetrials(trs - tr_rng(1) + 1,:,chs - ch_rng(1) + 1);
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

%summary of trial labels used
for k = 1:numel(tr_conditional_use)
    fprintf('Trial label "%s": %d trials flagged, %d file(s) not evaluated\n', ...
        tr_conditional_use{k}, trcond_nflagged(k), numel(trcond_noteval{k}))
    if ~isempty(trcond_noteval{k})
        warning('Trial label "%s" was never evaluated for (column is NaN for these): %s', ...
            tr_conditional_use{k}, strjoin(trcond_noteval{k}, ', '))
    end
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
    triallabel_trialnum,triallabel_region,triallabel_chanID,...
    'VariableNames', ["P1_MUAAUC","P2_MUAAUC","All_MUAAUC","P1_MUAPeak","P2_MUAPeak","P1_MUAPeakTime","P2_MUAPeakTime",...
    "Animal_Name","Animal_Num","Condition_Name","Condition_Num","Condition_FileNum","Trial","Region","Channel_ID"]);

%trial label columns, one per tr_conditional_use name
for k = 1:numel(tr_conditional_use)
    tableres.(matlab.lang.makeValidName("Excl_" + tr_conditional_use{k})) = triallabel_trcond(:,k);
end

%% add channel / shank label columns (from ProbeInfo; can be rerun on its own once tableres exists)
%Chan_Shank = shank (ProbeMaps column) of the channel
%Chan_Depth = position among the shank's kept (not Ch_Remove) channels, 1 = top (ProbeMaps row 1)
%ShankRank_<score> = rank of the channel's shank by that ShankLabels score, 1 = highest;
%   NaN if that mouse's probe has no such score
tableres.Chan_Shank = nan(height(tableres),1);
tableres.Chan_Depth = nan(height(tableres),1);
rankcols = matlab.lang.makeValidName("ShankRank_" + string(shankscore_use));
for k = 1:numel(rankcols)
    tableres.(rankcols(k)) = nan(height(tableres),1);
end

for m = unique(tableres.Animal_Num)'
    rows = find(tableres.Animal_Num == m);
    pinfo = load(fullfile(['E:\Roy\Processed Silicon Probe Data\ProbeInfo\' char(mice{m}) '-ProbeInfo.mat']),'ProbeInfo');
    prb = find(strcmp(pinfo.ProbeInfo.Areas, region{m}));
    assert(numel(prb) == 1, '%s: region "%s" matches %d probes in ProbeInfo.Areas (need exactly 1)', ...
        mice{m}, region{m}, numel(prb))
    pmap = pinfo.ProbeInfo.ProbeMaps{prb};
    chremove = [];
    if numel(pinfo.ProbeInfo.Ch_Remove) >= prb
        chremove = pinfo.ProbeInfo.Ch_Remove{prb};
    end

    %kept channels numbered 1..n from the top of each shank
    depthmap = nan(size(pmap));
    for s = 1:size(pmap,2)
        kept = pmap(:,s) > 0 & ~ismember(pmap(:,s), chremove);
        depthmap(kept,s) = 1:nnz(kept);
    end

    [onmap, loc] = ismember(tableres.Channel_ID(rows), pmap);
    assert(all(onmap), '%s: some Channel_IDs in the table are not on probe %d (%s)', mice{m}, prb, region{m})
    [~, shk] = ind2sub(size(pmap), loc);
    tableres.Chan_Shank(rows) = shk;
    tableres.Chan_Depth(rows) = depthmap(loc);
    if any(isnan(depthmap(loc)))
        warning('%s: some channels in the table are now in Ch_Remove; their Chan_Depth is NaN', mice{m})
    end

    %shank ranks: highest score = 1, NaN score -> NaN rank, ties -> lower shank number first
    SL = [];
    if isfield(pinfo.ProbeInfo,'ShankLabels') && numel(pinfo.ProbeInfo.ShankLabels) >= prb
        SL = pinfo.ProbeInfo.ShankLabels{prb};
    end
    for k = 1:numel(shankscore_use)
        sname = char(shankscore_use(k));
        if isstruct(SL) && isfield(SL, sname) && isnumeric(SL.(sname)) && numel(SL.(sname)) == size(pmap,2)
            score = double(SL.(sname));
            shankrank = nan(1,numel(score));
            valid = find(~isnan(score));
            [~,ord] = sort(score(valid),'descend');
            shankrank(valid(ord)) = 1:numel(valid);
            tableres.(rankcols(k))(rows) = shankrank(shk);
        end
    end
end
clear pinfo SL

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
