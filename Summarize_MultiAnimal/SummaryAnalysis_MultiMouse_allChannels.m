% SummaryAnalysis_MultiMouse_allChannels.m
%
% Description: Combined all-channels summary of LFP, MUA and time-frequency (TF)
%   responses across mice and conditions, in one long-format table (one row per
%   good trial x kept channel, all measures side by side) written to a csv for R.
%   Merges SummaryAnalysis_MultiMouseLFP_allChannels.m,
%   SummaryAnalysis_MultiMouseMUA_modifiedforallChannels.m and
%   SummaryAnalysis_MultiMouseTF_allChannels.m (which remain available for one data
%   type at a time): same start and end sections, one calculation block per data
%   type in the middle. All three types are always computed.
%   The user lists the mice, the region per mouse, the condition names and
%   (optionally) trial labels, picks shank activity scores and then only the LFP
%   files. Each recording's spiking file is found from the LFP file name; TF is
%   recomputed from the LFP (no TF files are read). The probe is the one whose
%   ProbeInfo.Areas matches the mouse's region; its channels minus
%   ProbeInfo.Ch_Remove are analyzed.
%   LFP: area under the RMS envelope of the lfp_band band (zero-phase Butterworth,
%     each trial filtered on its own), and the negative peak and its latency on a
%     peak_smooth_ms moving mean (no threshold), in P1 and P2 windows.
%   MUA: spikes binned to 1 ms and smoothed into a firing rate (50 ms moving mean);
%     AUC in P1/P2/All windows, and the peak rate and its latency in P1/P2.
%   TF: Morlet wavelet power in each band (whole-Hz frequencies, the same cycle rule
%     as OpenEphys_BaseAnalysis*.m), in P1/P2/All windows relative to the same
%     trial's baseline.
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
%     lfp_band (double, 1 x 2) - bandpass for the LFP RMS envelope, Hz ([1 150]); Butterworth,
%       2nd order per edge, second-order sections, filtfilt (zero phase)
%     rms_win_ms (double) - RMS envelope window for LFP RMSsum (25)
%     peak_smooth_ms (double) - moving mean used to find LFP peaks (11)
%     fs_mua (double) - sampling rate of the spiking data, Hz (30000)
%     tf_bands (cell, nBands x 2) - band name (string, used in the column names) and the list of
%       frequencies in Hz that make up the band (e.g. 10:18)
%   Files picked per condition via uipickfiles (superCondis_dir, cell 1 x nCond), from
%   E:\Roy\Processed Silicon Probe Data\LFP\<animal>\<animal>-<stim>-LFP.mat, each containing:
%     stim_lfp_stimchunks (double, trials x samples x channels) - peri-stimulus LFP, 1 kHz,
%       presumably in uV [inferred: units not stated upstream]. Stim onset is sample 5001
%       (stim_onset): OpenEphys_BaseAnalysis*.m epochs trials as stim_times-pre:stim_times+post
%       with pre = 5000 samples. Channel index is assumed to equal the raw channel ID
%       [inferred: holds when ProbeInfo.ChanIds = 1:Chans]
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1 of the data)
%     tr_remove (0/1, 1 x nTrials) - mask over tr_keep (updated format, as written by
%       ConditionalTrialRemove_callable.m); only tr_keep(~tr_remove) are analyzed.
%       Must be the same length as tr_keep, else the script stops with an error
%     tr_remove_conditional (table, columns Name (string) / TrialIdx (double), optional) -
%       trial labels from ConditionalTrialRemove_callable.m; one Name can have many rows;
%       NaN TrialIdx = evaluated, none flagged; Name absent = file not evaluated for it
%   Spiking file of each recording, found by replacing \LFP\ with \Spiking\ and -LFP.mat with
%   -spiking_results.mat in the LFP file's path, containing:
%     stim_spike_stimchunks (0/1, trials x samples x channels) - spike indicator at 30 kHz,
%       same epochs as the LFP (stim onset in 1 ms bin 5001)
%     tr_keep, tr_remove - must equal the LFP file's (the script stops otherwise)
%     If a spiking file is missing, a dialog lists them and offers abort or continue
%     (continue: those recordings get NaN in the MUA columns)
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
%   spikefile, hasMUA (cell, 1 x nCond) - per file, spiking file path and whether it exists
%   mouseProbe (double, 1 x nMice) - probe index used for each mouse
%   chanlist (cell, 1 x nMice) - sorted column vector of raw channel IDs analyzed per mouse
%   stim_onset, P1wind, P2wind, Allwind (double) - LFP/MUA stimulus sample (5001) and windows
%     as samples: 75-350, 350-3000 and 75-3000 ms post stimulus
%   tf_seg, tf_onset, tfP1wind, tfP2wind, tfAllwind, tfbasewind (double) - TF segment (LFP
%     samples 2000:8000, as OpenEphys_BaseAnalysis*.m), stimulus sample within it (3002),
%     windows 75-350, 350-2999, 75-2999 ms and baseline -2102..-102 ms
%   lfp_sos, lfp_g (double) - the lfp_band Butterworth filter as second-order sections and gain
%   metriccols (struct) - one field per measure (column name), each a column over all rows
%
% Outputs:
%   tableres (table, nRows x (13 + 3*nBands + 13 + nLabels + nScores)) - one row per good trial x
%   kept channel, ordered condition -> file -> channel -> trial. Columns:
%     LFP (measured on the LFP file):
%     P1_LFPRMSsum, P2_LFPRMSsum (double) - area under the rms_win_ms RMS envelope of the
%       lfp_band-filtered LFP over P1 (75-350 ms post stimulus) / P2 (350-3000 ms), in uV*ms.
%       The filter smears ~0.3 s, so a strong P1 response adds ~12% of its own RMSsum to P2
%     P1_LFPPeak, P2_LFPPeak (double) - size of the negative peak in the P1 / P2 window, as a
%       positive value in uV, measured on a peak_smooth_ms moving mean of the unfiltered LFP:
%       the highest local maximum of the sign-flipped trace (lower on both sides, window padded
%       1 sample each side; a flat top counts once, at its center). No threshold: noise-only
%       trials also get a (small) peak. NaN only if the window has no local max
%     P1_LFPPeakTime, P2_LFPPeakTime (double) - time of that peak, ms post stimulus; NaN when Peak is NaN
%     MUA (measured on the spiking file; NaN if it is missing):
%     P1_MUAAUC, P2_MUAAUC, All_MUAAUC (double) - integral of the firing rate over P1, P2 and
%       All (75-3000 ms), in spikes/s*ms
%     P1_MUAPeak, P2_MUAPeak (double) - rate at the highest local maximum in the P1 / P2 window
%       (lower on both sides, window padded 1 bin each side; a flat top counts once, at its
%       center), in spikes/s. 0 if the rate is 0 throughout the window; NaN if there are spikes
%       but no local max (rate only decaying/rising)
%     P1_MUAPeakTime, P2_MUAPeakTime (double) - time of that peak, ms post stimulus; NaN
%       whenever there is no peak (Peak 0 or NaN)
%     TF (recomputed from the LFP file):
%     P1_TF<band>_Pow, P2_TF<band>_Pow, All_TF<band>_Pow (double, 3 per tf_bands row) - mean
%       Morlet wavelet power over the band's frequencies and the P1 (75-350 ms), P2 (350-2999 ms)
%       or All (75-2999 ms) window, divided by the same trial's mean band power over the
%       baseline (-2102..-102 ms). Ratio, unitless; 1 = no change. Wavelet as in
%       OpenEphys_BaseAnalysis*.m (4 s long, all trials of the file joined end to end, LFP
%       samples 2000:8000, single precision) with cycles(f) = 6*(10/6)^((f-1)/149), the same
%       6-to-10-cycle rule BaseAnalysis uses over 1-150 Hz; on whole-Hz frequencies, so values
%       differ very slightly from the stored stim_tf (which uses 0.993 Hz steps)
%     Animal_Name (string), Animal_Num (double) - mouse ID and its index in mice
%     Condition_Name (string), Condition_Num (double) - condition name and its index in condis
%     Condition_FileNum (double) - Nth file of this mouse within this condition (1..n)
%     Trial (double) - original trial number (value from tr_keep)
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
% Dependencies: OpenEphys_BaseAnalysis*.m (writes *-LFP.mat, *-spiking_results.mat and
%   <animal>-ProbeInfo.mat); OpenEphys_editProbeInfo_ChRemove.m (sets Ch_Remove);
%   ConditionalTrialRemove.m / ConditionalTrialRemove_callable.m (tr_remove mask format,
%   tr_remove_conditional labels); ProbeInfo.ShankLabels from the channel/shank labelling
%   scripts. Requires uipickfiles (File Exchange), the Signal Processing Toolbox (butter,
%   zp2sos, filtfilt, envelope) and MATLAB R2017b+ (islocalmax).

%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%names must match the start of the data file names (animal-stim-LFP.mat)
mice = {"20260423-p12"};%, "20260402-p12"};

%which brain area is being recorded for each mouse? (same order as mice)
%used to pick the probe: must match that mouse's ProbeInfo.Areas
region = {"V1"};

condis = {'W','L_4','LW_4', 'L_8','LW_8','L_12','LW_12','L_15','LW_15'}; %list conditions to be included, named as you'd prefer. Note the order
% as they will be treated as factors later

%trial labels made with ConditionalTrialRemove.m to add as columns (empty = none)
%each becomes a column Excl_<name>: 1 = flagged, 0 = evaluated & not flagged, NaN = file not evaluated
tr_conditional_use = {}; %e.g. {"whisker_twitch_artifact"}

%data types to include - not used yet: all three are always computed
%(the single-type scripts SummaryAnalysis_MultiMouseLFP/MUA/TF_allChannels.m cover other cases)
%do_LFP = true;
%do_MUA = true;
%do_TF = true;

%LFP settings
lfp_band = [1 150];   %Hz, bandpass for the RMS envelope
rms_win_ms = 25;      %RMS envelope window for RMSsum
peak_smooth_ms = 11;  %moving mean used to find peaks

%MUA settings
fs_mua = 30000;       %Hz, sampling rate of stim_spike_stimchunks

%TF settings
tf_bands = {"AlphaBeta",10:18;   %band name (goes into the column names), frequencies in Hz
            "LoloGamma",30:50};  %add rows for more bands

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
%% pick the LFP files for each condition (spiking files are found from their names)
for con = 1:length(condis)
    superCondis_dir{con} = uipickfiles('FilterSpec','E:\Roy\Processed Silicon Probe Data\LFP','Prompt', ['Choose ' condis{con} ' LFP files']);
end
%% match files to mice, find each recording's spiking file, and get channels for each mouse from ProbeInfo

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

%spiking file of each recording (for MUA): same folder structure and name as the LFP file
spikefile = cell(1,numel(condis));
hasMUA = cell(1,numel(condis));
missingspk = strings(0,1);
for con = 1:numel(condis)
    for file = 1:numel(superCondis_dir{con})
        lf = superCondis_dir{con}{file};
        assert(endsWith(lf,'-LFP.mat'), 'Condition %d (%s): "%s" is not an LFP file (<animal>-<stim>-LFP.mat)', con, condis{con}, lf)
        spikefile{con}{file} = strrep(strrep(lf,'\LFP\','\Spiking\'),'-LFP.mat','-spiking_results.mat');
        hasMUA{con}(file) = isfile(spikefile{con}{file});
        if ~hasMUA{con}(file)
            missingspk(end+1,1) = string(spikefile{con}{file}); %#ok<SAGROW>
        end
    end
end
if ~isempty(missingspk)
    warning('No spiking file for %d recording(s):\n%s', numel(missingspk), strjoin(missingspk, newline))
    answer = questdlg(sprintf('No spiking file found for %d recording(s) (their MUA columns will be NaN):\n\n%s', ...
        numel(missingspk), strjoin(missingspk, newline)), 'Missing spiking files', 'Continue (MUA = NaN)', 'Abort', 'Abort');
    if ~strcmp(answer,'Continue (MUA = NaN)')
        error('Aborted by user: spiking files are missing for some recordings')
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
%% load each recording and compute the LFP, MUA and TF measures, one file at a time
metriccols = struct(); %one field per measure (column name), filled by the calculation blocks

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

%LFP / MUA timing: 1 ms samples (LFP) and 1 ms bins (MUA)
stim_onset = 5001; %sample of stimulus onset (0 ms post stimulus); epochs are stim_times-pre:stim_times+post
P1wind = stim_onset + (75:350);   %ms post stimulus 75-350
P2wind = stim_onset + (350:3000); %350-3000
Allwind = stim_onset + (75:3000); %75-3000

%LFP bandpass for the RMS envelope: Butterworth, 2nd order per edge, as second-order sections (stable
%at low cutoffs), applied forward and backward with filtfilt (zero phase). LFP is 1 kHz
[z,p,k] = butter(2, lfp_band/(1000/2), 'bandpass');
[lfp_sos,lfp_g] = zp2sos(z,p,k);

%TF: the wavelet is run on LFP samples 2000:8000 of all trials joined end to end, as in OpenEphys_BaseAnalysis*.m
tf_seg = 2000:8000;
tf_onset = stim_onset - tf_seg(1) + 1;   %stimulus sample within the TF segment (3002)
tfP1wind = tf_onset + (75:350);          %ms post stimulus 75-350
tfP2wind = tf_onset + (350:2999);        %350-2999 (last sample of the segment is +2999 ms)
tfAllwind = tf_onset + (75:2999);        %75-2999
tfbasewind = tf_onset + (-2102:-102);    %baseline, as in the TF scripts
wavtime = -2:1/1000:2;                   %wavelet time axis, s (4 s long, as BaseAnalysis)
halfwave = (numel(wavtime)-1)/2;
tf_cycles = @(f) 6*(10/6).^((f-1)/149);  %BaseAnalysis rule: 6 -> 10 cycles over 1-150 Hz (logarithmic)

%MUA columns (filled with NaN when a recording has no spiking file)
muacols = {'P1_MUAAUC','P2_MUAAUC','All_MUAAUC','P1_MUAPeak','P1_MUAPeakTime','P2_MUAPeak','P2_MUAPeakTime'};

for con = 1:numel(condis)
    for file = 1:numel(superCondis_dir{con})
        disp(['Running Condition ', num2str(con), ' (', condis{con}, '), File ', num2str(file)])
        m = condiFileMice{con}(file);
        chs = chanlist{m};
        dat = matfile(superCondis_dir{con}{file});
        [~,fname] = fileparts(superCondis_dir{con}{file});
        filevars = who(dat);

        %good trials (from the LFP file): tr_remove is a 0/1 mask over tr_keep
        tr_keep = dat.tr_keep;
        tr_remove = dat.tr_remove;
        %TEMPORARY: convert the old format (tr_keep = kept trial numbers, tr_remove = removed trial numbers
        %or empty) to the mask format (tr_keep = all trials, tr_remove = 0/1 mask). Remove once all files
        %have been through ConditionalTrialRemove(_callable).m
        if numel(tr_remove) ~= numel(tr_keep)
            alltr = sort([tr_keep(:); tr_remove(:)])';
            tr_remove = double(ismember(alltr, tr_remove));
            tr_keep = alltr;
            fprintf('%s: old tr_keep/tr_remove format converted to a mask (%d of %d trials removed)\n', fname, sum(tr_remove), numel(tr_keep))
        end
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

        % ============ LOAD ============
        %LFP: all trials (the TF wavelet sees them joined end to end), kept channels. LFP files are not
        %saved -v7.3, so the whole variable is loaded either way
        L = load(superCondis_dir{con}{file},'stim_lfp_stimchunks');
        assert(size(L.stim_lfp_stimchunks,1) == numel(tr_keep), '%s: stim_lfp_stimchunks has %d trials but tr_keep has %d', ...
            fname, size(L.stim_lfp_stimchunks,1), numel(tr_keep))
        lfpall = L.stim_lfp_stimchunks(:,:,chs);
        clear L
        lfpgood = lfpall(trs,:,:); %good trials x samples x kept channels

        %MUA: spiking file of the same recording; same trial mask required.
        %matfile can only read evenly spaced ranges, so read the block spanning the good trials
        %and kept channels once, then pick the ones needed from it in memory
        if hasMUA{con}(file)
            sdat = matfile(spikefile{con}{file});
            sk = sdat.tr_keep;  sr = sdat.tr_remove;
            if numel(sr) ~= numel(sk) %TEMPORARY: same old-format conversion as for the LFP file
                allsk = sort([sk(:); sr(:)])';
                sr = double(ismember(allsk, sr));
                sk = allsk;
            end
            assert(isequal(sk(:),tr_keep(:)) && isequal(sr(:),tr_remove(:)), ...
                '%s: tr_keep/tr_remove in the spiking file differ from the LFP file', fname)
            tr_rng = min(trs):max(trs);
            ch_rng = min(chs):max(chs);
            spkdata = sdat.stim_spike_stimchunks(tr_rng,:,ch_rng);
            spkdata = spkdata(trs - tr_rng(1) + 1,:,chs - ch_rng(1) + 1);
            nbatch = floor(size(spkdata,2)/(fs_mua/1000));
        end

        %TF: Morlet wavelet power of each band, for the kept channels and good trials, computed as
        %OpenEphys_BaseAnalysis*.m builds stim_tf (all trials joined, single precision)
        tfdata = single(lfpall(:,tf_seg,:));
        ntrall = size(tfdata,1);
        ntfsamp = numel(tf_seg);
        nConv = numel(wavtime) + ntfsamp*ntrall - 1;
        dX = fft(reshape(permute(tfdata,[2 1 3]),ntfsamp*ntrall,nch),nConv); %columns = channels, trials joined
        clear tfdata lfpall
        bandpow = cell(1,size(tf_bands,1)); %each: channels x samples x good trials
        for b = 1:size(tf_bands,1)
            acc = zeros(nch,ntfsamp,ntr);
            for fq = tf_bands{b,2}
                wwidth = tf_cycles(fq)/(2*pi*fq);
                wX = fft(exp(2*1i*pi*fq.*wavtime) .* exp(-wavtime.^2./(2*wwidth^2)),nConv).';
                as = ifft((wX./max(wX)) .* dX);
                as = reshape(as(halfwave+1:end-halfwave,:),ntfsamp,ntrall,nch);
                acc = acc + double(permute(abs(as(:,trs,:)).^2,[3 1 2]));
            end
            bandpow{b} = acc / numel(tf_bands{b,2});
        end
        clear dX as wX acc
        % ==============================

        for ch = 1:nch
            chmetrics = struct(); %one field per measure, each ntr x 1 (one value per trial)

            % ============ LFP ============
            chdata = reshape(lfpgood(:,:,ch),ntr,[]); %good trials x samples for this channel
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
            % =============================

            % ============ MUA ============
            if hasMUA{con}(file)
                chspikes = reshape(spkdata(:,:,ch),ntr,[]);
                spikeper = zeros(ntr,nbatch);
                for batch = 1:nbatch
                    spikebatchi = sum(chspikes(:,1+((fs_mua/1000)*(batch-1)):(fs_mua/1000)+((fs_mua/1000)*(batch-1))),2)>0;
                    spikeper(:,batch) = spikebatchi;
                end
                spikes_rate = movmean(spikeper,50,2)*1000; %spikes/s, good trials x 1 ms bins

                auc = cumtrapz(spikes_rate(:,P1wind),2);
                chmetrics.P1_MUAAUC = auc(:,end);
                auc = cumtrapz(spikes_rate(:,P2wind),2);
                chmetrics.P2_MUAAUC = auc(:,end);
                auc = cumtrapz(spikes_rate(:,Allwind),2);
                chmetrics.All_MUAAUC = auc(:,end);

                %Peaks: highest local max in the window (lower on both sides), and its time in ms post stimulus.
                %window is padded 1 bin each side so a decay into / rise out of the window isn't a peak.
                %rounded so movmean plateaus are exactly flat, and a plateau's peak is its center.
                %no local max (only decay/rise) -> Peak and PeakTime NaN; no spikes (rate all 0) -> Peak 0, PeakTime NaN
                seg = round(spikes_rate(:,P1wind(1)-1:P1wind(end)+1),6);
                lm = islocalmax(seg,2,'FlatSelection','center');
                lm(:,[1 end]) = false;
                seg(~lm) = -Inf;
                [peaks,tps] = max(seg,[],2);
                peaktime = tps + P1wind(1) - 2 - stim_onset; %padded index -> bin -> ms post stimulus
                nopeak = isinf(peaks);
                allzero = all(spikes_rate(:,P1wind) == 0,2);
                peaks(nopeak) = NaN;  peaks(allzero) = 0;
                peaktime(nopeak | allzero) = NaN;
                chmetrics.P1_MUAPeak = peaks;
                chmetrics.P1_MUAPeakTime = peaktime;

                seg = round(spikes_rate(:,P2wind(1)-1:P2wind(end)+1),6);
                lm = islocalmax(seg,2,'FlatSelection','center');
                lm(:,[1 end]) = false;
                seg(~lm) = -Inf;
                [peaks,tps] = max(seg,[],2);
                peaktime = tps + P2wind(1) - 2 - stim_onset; %padded index -> bin -> ms post stimulus
                nopeak = isinf(peaks);
                allzero = all(spikes_rate(:,P2wind) == 0,2);
                peaks(nopeak) = NaN;  peaks(allzero) = 0;
                peaktime(nopeak | allzero) = NaN;
                chmetrics.P2_MUAPeak = peaks;
                chmetrics.P2_MUAPeakTime = peaktime;
            else
                for q = 1:numel(muacols)
                    chmetrics.(muacols{q}) = nan(ntr,1); %no spiking file for this recording
                end
            end
            % =============================

            % ============ TF ============
            %band power in each window / the same trial's baseline band power (1 = no change)
            for b = 1:size(tf_bands,1)
                pw = reshape(permute(bandpow{b}(ch,:,:),[3 2 1]),ntr,[]); %trials x samples
                base = mean(pw(:,tfbasewind),2);
                chmetrics.(char("P1_TF" + tf_bands{b,1} + "_Pow")) = mean(pw(:,tfP1wind),2) ./ base;
                chmetrics.(char("P2_TF" + tf_bands{b,1} + "_Pow")) = mean(pw(:,tfP2wind),2) ./ base;
                chmetrics.(char("All_TF" + tf_bands{b,1} + "_Pow")) = mean(pw(:,tfAllwind),2) ./ base;
            end
            % ============================

            %append this channel's measures to metriccols (any column names; same set every channel)
            fns = fieldnames(chmetrics);
            if ~isempty(fieldnames(metriccols))
                assert(isempty(setxor(fns, fieldnames(metriccols))), ...
                    '%s: the calculation blocks must set the same measures for every channel', fname)
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
        clear lfpgood spkdata bandpow chdata chenv chsm pw
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
