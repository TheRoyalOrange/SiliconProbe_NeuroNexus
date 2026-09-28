% SummaryAnalysis_MultiMouseLFP_allChannels.m
%
% Description: Pools peri-stimulus LFP across mice and conditions using every
%   usable channel on one probe per mouse (same structure as
%   SummaryAnalysis_MultiMouseMUA_modifiedforallChannels.m). The user lists only
%   the mice, the region per mouse, the condition names and (optionally) trial
%   labels, picks shank activity scores and then the files. The probe is the one
%   whose ProbeInfo.Areas matches the mouse's region; its channels minus
%   ProbeInfo.Ch_Remove are analyzed. Each file's mouse is found from its file
%   name. For every good trial x channel, the measures set in the
%   DATA-TYPE-SPECIFIC CALCULATION block are computed, then all rows go into one
%   long-format table (with trial, channel, shank/depth and label columns) that
%   is written to a csv for R.
%   LFP measures (ported from SummaryAnalysis_MultiMouseLFP.m): area under the
%   RMS envelope of the lfp_band (1-150 Hz) band, and the negative peak and its
%   latency, in an early (P1) and late (P2) window after the stimulus. Each trial
%   is filtered on its own with a zero-phase Butterworth bandpass (the old script's
%   steep IIR bandpass distorted the signal - a 60 uV dip came out ~2x larger with
%   ringing across the whole pre-stim period - so RMSsum is NOT comparable to CSVs
%   from the old script). Peaks are found on a peak_smooth_ms moving mean of the
%   unfiltered LFP, with no thresholding: every window with a local max gets a peak.
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
%     datadir (char) - folder the result files are picked from
%     datavar (char) - name of the data variable in each result file
%     lfp_band (double, 1 x 2) - bandpass for the RMS envelope, Hz ([1 150]); Butterworth,
%       2nd order per edge, second-order sections, filtfilt (zero phase)
%     rms_win_ms (double) - RMS envelope window for RMSsum (25)
%     peak_smooth_ms (double) - moving mean used to find peaks (11)
%   Files picked per condition via uipickfiles (superCondis_dir, cell 1 x nCond),
%   from <datadir>\<animal>\<animal>-<stim>-LFP.mat, each containing:
%     stim_lfp_stimchunks (double, trials x samples x channels) - peri-stimulus LFP,
%       1 kHz, presumably in uV [inferred: units not stated upstream]. Stim onset is
%       sample 5001 (stim_onset): OpenEphys_BaseAnalysis*.m epochs trials as
%       stim_times-pre:stim_times+post with pre = 5000 samples, and stim_times is the
%       first high TTL sample. Channel index is assumed to equal the raw channel ID
%       [inferred: holds when ProbeInfo.ChanIds = 1:Chans]
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1 of the data)
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
%     ProbeInfo.Ch_Remove (cell, 1 x nProbes) - raw channel IDs to exclude per probe;
%       missing/empty = exclude none. Assumed indexed by probe number
%       [inferred: OpenEphys_BaseAnalysis*.m indexes it by position in poi, which is
%       the same only when poi = 1:ProbeNum (the default)]
%     ProbeInfo.ShankLabels (cell, 1 x nProbes, optional) - per probe a struct (or [] if not labelled);
%       each field holding one number per shank (1 x nShanks, ProbeMaps column order, raw values)
%       is offered as a shank activity score (e.g. SponActivity, LightOnly_mixedIntensity)
%
% Derived (not user-set):
%   shankscore_use (string, 1 x nScores) - shank activity scores chosen in a list dialog at the
%     start (from those found across the mice); if some mice lack one, a dialog offers
%     abort or continue (those mice get NaN in that column)
%   condiFileMice (cell, 1 x nCond) - per file, index into mice (from the file name)
%   condiFileNum (cell, 1 x nCond) - per file, Nth file of that mouse in that condition (1..n)
%   mouseProbe (double, 1 x nMice) - probe index used for each mouse
%   chanlist (cell, 1 x nMice) - sorted column vector of raw channel IDs analyzed per mouse
%   stim_onset, P1wind, P2wind, Allwind (double) - stimulus sample (5001) and the response
%     windows as samples: 75-350, 350-3000 and 75-3000 ms post stimulus
%   metriccols (struct) - one field per measure (column name), each a column over all rows
%   lfp_sos, lfp_g (double) - the lfp_band Butterworth filter as second-order sections and gain
%
% Outputs:
%   tableres (table, nRows x (nMeasures + 11 + nLabels + nScores)) - one row per good
%   trial x kept channel, ordered condition -> file -> channel -> trial. Columns:
%     (measure columns, set in the calculation block; names carry "LFP" so LFP/MUA/TF
%     tables can be combined:)
%     P1_LFPRMSsum, P2_LFPRMSsum (double) - area under the rms_win_ms RMS envelope of the
%       lfp_band-filtered LFP over P1 (75-350 ms post stimulus) / P2 (350-3000 ms), in uV*ms.
%       The filter smears ~0.3 s, so a strong P1 response adds ~12% of its own RMSsum to P2
%     P1_LFPPeak, P2_LFPPeak (double) - size of the negative peak in the P1 / P2 window, as a
%       positive value in uV, measured on a peak_smooth_ms moving mean of the unfiltered LFP:
%       the highest local maximum of the sign-flipped trace (lower on both sides, window padded
%       1 sample each side so a deflection still decaying into / rising out of the window isn't
%       a peak; a flat top counts once, at its center). No threshold: noise-only trials also get
%       a (small) peak. NaN only if the window has no local max
%     P1_LFPPeakTime, P2_LFPPeakTime (double) - time of that peak, in ms post stimulus
%       (sample - stim_onset); NaN when Peak is NaN
%     Animal_Name (string), Animal_Num (double) - mouse ID and its index in mice
%     Condition_Name (string), Condition_Num (double) - condition name and its index in condis
%     Condition_FileNum (double) - Nth file of this mouse within this condition (1..n)
%     Trial (double) - original trial number (value from tr_keep); with Animal/Condition/
%       FileNum/Channel_ID, identifies a row across the LFP/MUA/TF tables
%     Region (string) - from region
%     Channel_ID (double) - raw channel ID (unique across probes)
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
% Dependencies: OpenEphys_BaseAnalysis*.m (writes *-LFP.mat and <animal>-ProbeInfo.mat);
%   OpenEphys_editProbeInfo_ChRemove.m (sets Ch_Remove);
%   ConditionalTrialRemove.m / ConditionalTrialRemove_callable.m (tr_remove mask format,
%   tr_remove_conditional labels); ProbeInfo.ShankLabels from the channel/shank labelling
%   scripts. Requires uipickfiles (File Exchange), the Signal Processing Toolbox (butter,
%   zp2sos, filtfilt, envelope) and MATLAB R2017b+ (islocalmax).

%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%names must match the start of the data file names (animal-stim-LFP.mat)
mice = {"20260226-p12"};%, "20260402-p12"};

%which brain area is being recorded for each mouse? (same order as mice)
%used to pick the probe: must match that mouse's ProbeInfo.Areas
region = {"V1"};

condis = {'L_4','LW_4', 'L_8','LW_8','L_12','LW_12','L_15','LW_15'}; %list conditions to be included, named as you'd prefer. Note the order
% as they will be treated as factors later

%trial labels made with ConditionalTrialRemove.m to add as columns (empty = none)
%each becomes a column Excl_<name>: 1 = flagged, 0 = evaluated & not flagged, NaN = file not evaluated
tr_conditional_use = {}; %e.g. {"whisker_twitch_artifact"}

%data type (the only place it is named)
datadir = 'E:\Roy\Processed Silicon Probe Data\LFP'; %where the result files are picked from
datavar = 'stim_lfp_stimchunks';                      %data variable in each file: trials x samples x channels
lfp_band = [1 150];   %Hz, bandpass for the RMS envelope
rms_win_ms = 25;      %RMS envelope window for RMSsum
peak_smooth_ms = 11;  %moving mean used to find peaks

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
%% pick the files for each condition
for con = 1:length(condis)
    superCondis_dir{con} = uipickfiles('FilterSpec',datadir,'Prompt', ['Choose ' condis{con} ' files']);
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
%% load each file and compute the measures, one file at a time
metriccols = struct(); %one field per measure (column name), filled by the calculation block

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

stim_onset = 5001; %sample (1 ms) of stimulus onset (0 ms post stimulus); epochs are stim_times-pre:stim_times+post
P1wind = stim_onset + (75:350);   %ms post stimulus 75-350
P2wind = stim_onset + (350:3000); %350-3000
Allwind = stim_onset + (75:3000); %75-3000

%bandpass for the RMS envelope: Butterworth, 2nd order per edge, as second-order sections (stable at
%low cutoffs), applied forward and backward with filtfilt (zero phase). LFP is 1 kHz
[z,p,k] = butter(2, lfp_band/(1000/2), 'bandpass');
[lfp_sos,lfp_g] = zp2sos(z,p,k);

for con = 1:numel(condis)
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
        filedata = dat.(datavar)(tr_rng,:,ch_rng);
        filedata = filedata(trs - tr_rng(1) + 1,:,chs - ch_rng(1) + 1);

        for ch = 1:nch
            chdata = reshape(filedata(:,:,ch),ntr,[]); %trials x samples for this channel
            chmetrics = struct(); %one field per measure, each ntr x 1 (one value per trial)

            % ============ DATA-TYPE-SPECIFIC CALCULATION ============
            %compute the measures from chdata (good trials x samples, this channel) and set
            %chmetrics.<ColumnName> = ntr x 1 column for each, e.g. chmetrics.P1_LFPRMSsum
            %(windows: P1wind / P2wind / Allwind; latency in ms = sample - stim_onset)

            %RMS envelope of the lfp_band band; each trial filtered on its own, zero phase (columns of chdata')
            chenv = envelope(filtfilt(lfp_sos,lfp_g,chdata'),rms_win_ms,'rms')';
            %sign-flipped moving mean used to find negative peaks (negative deflections become maxima)
            chsm = -movmean(chdata,peak_smooth_ms,2);

            %For each window: RMSsum (area under the envelope), and the negative peak = highest local
            %max of chsm in the window (lower on both sides; window padded 1 sample each side so a
            %deflection still decaying into / rising out of it isn't a peak). Peak size is from the
            %smoothed trace, as a positive value; NaN (with its time) if the window has no local max.

            %P1
            auc = cumtrapz(chenv(:,P1wind),2);
            chmetrics.P1_LFPRMSsum = auc(:,end);
            seg = chsm(:,P1wind(1)-1:P1wind(end)+1);
            lm = islocalmax(seg,2,'FlatSelection','center');
            lm(:,[1 end]) = false;
            seg(~lm) = -Inf;
            [peaks,tps] = max(seg,[],2);
            peaktime = tps + P1wind(1) - 2 - stim_onset; %padded index -> sample -> ms post stimulus
            nopeak = isinf(peaks);
            peaks(nopeak) = NaN;  peaktime(nopeak) = NaN;
            chmetrics.P1_LFPPeak = peaks;
            chmetrics.P1_LFPPeakTime = peaktime;

            %P2
            auc = cumtrapz(chenv(:,P2wind),2);
            chmetrics.P2_LFPRMSsum = auc(:,end);
            seg = chsm(:,P2wind(1)-1:P2wind(end)+1);
            lm = islocalmax(seg,2,'FlatSelection','center');
            lm(:,[1 end]) = false;
            seg(~lm) = -Inf;
            [peaks,tps] = max(seg,[],2);
            peaktime = tps + P2wind(1) - 2 - stim_onset; %padded index -> sample -> ms post stimulus
            nopeak = isinf(peaks);
            peaks(nopeak) = NaN;  peaktime(nopeak) = NaN;
            chmetrics.P2_LFPPeak = peaks;
            chmetrics.P2_LFPPeakTime = peaktime;
            % =========================================================

            %append this channel's measures to metriccols (any column names; same set every channel)
            fns = fieldnames(chmetrics);
            if ~isempty(fieldnames(metriccols))
                assert(isempty(setxor(fns, fieldnames(metriccols))), ...
                    '%s: the calculation block must set the same measures for every channel', fname)
            end
            for f = 1:numel(fns)
                assert(isequal(size(chmetrics.(fns{f})),[ntr 1]), ...
                    '%s: measure %s must be ntr x 1 (one value per trial)', fname, fns{f})
                if ~isfield(metriccols,fns{f})
                    metriccols.(fns{f}) = [];
                end
                metriccols.(fns{f}) = [metriccols.(fns{f}); chmetrics.(fns{f})];
            end
        end
        clear filedata chdata chenv chsm
    end
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
%% make a table that can be easily read out in R for stats and plotting
labeltab = table(triallabel_animalname, triallabel_animalnum, triallabel_condiname,triallabel_condinum,triallabel_condifile,...
    triallabel_trialnum,triallabel_region,triallabel_chanID,...
    'VariableNames', ["Animal_Name","Animal_Num","Condition_Name","Condition_Num","Condition_FileNum","Trial","Region","Channel_ID"]);

%measure columns first (from metriccols), then the label columns
tableres = labeltab(:,[]);
fns = fieldnames(metriccols);
for f = 1:numel(fns)
    assert(numel(metriccols.(fns{f})) == height(labeltab), ...
        'Measure %s has %d rows but there are %d label rows', fns{f}, numel(metriccols.(fns{f})), height(labeltab))
    tableres.(fns{f}) = metriccols.(fns{f});
end
tableres = [tableres, labeltab];

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
