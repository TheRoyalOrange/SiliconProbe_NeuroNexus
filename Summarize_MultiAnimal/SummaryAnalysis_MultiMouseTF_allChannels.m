% SummaryAnalysis_MultiMouseTF_allChannels.m
%
% Description: Pools peri-stimulus time-frequency (Morlet wavelet) power across mice
%   and conditions using every usable channel on one probe per mouse (same
%   structure as SummaryAnalysis_MultiMouseMUA_modifiedforallChannels.m and
%   SummaryAnalysis_MultiMouseLFP_allChannels.m). The user lists only
%   the mice, the region per mouse, the condition names and (optionally) trial
%   labels, picks shank activity scores and then the files. The probe is the one
%   whose ProbeInfo.Areas matches the mouse's region; its channels minus
%   ProbeInfo.Ch_Remove are analyzed. Each file's mouse is found from its file
%   name. For every good trial x channel, the measures set in the
%   DATA-TYPE-SPECIFIC CALCULATION block are computed, then all rows go into one
%   long-format table (with trial, channel, shank/depth and label columns) that
%   is written to a csv for R.
%   TF measures (adapted from SummaryAnalysis_MultiMouseTF.m): for each band in
%   tf_bands, the mean band power in an early (P1), late (P2) and whole (All)
%   window after the stimulus, divided by the mean band power in that trial's
%   pre-stimulus baseline (a ratio; 1 = no change).
%   Speed: by default (tf_recompute) the band rows are recomputed from the same
%   recording's LFP file with exactly the method OpenEphys_BaseAnalysis*.m uses to
%   build stim_tf (matches the stored values to ~1e-6), which takes seconds per file.
%   Each file is checked against 3 stored stim_tf samples first; if they differ (e.g.
%   custom TF settings were used) or the LFP file is missing, stim_tf is read instead.
%   Reading is slow: the TF files are 15-42 GB and stored in chunks that each hold all
%   channels and trials for one frequency at one sample, so reading one channel costs
%   as much as all of them (~12 s per frequency row, ~5-6 min per file for the default
%   bands). The reading path reads all kept channels at once, one frequency row at a
%   time, turning each into power straight away (low memory).
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
%     tf_bands (cell, nBands x 2) - band name (string, used in the column names) and
%       [low high] frequency range in Hz, inclusive; add rows for more bands
%     tf_recompute (logical) - true: recompute the band rows from the LFP file (fast, checked
%       against stim_tf per file, falls back to reading it); false: always read stim_tf
%   Files picked per condition via uipickfiles (superCondis_dir, cell 1 x nCond),
%   from <datadir>\<animal>\<animal>-<stim>_TF_results.mat, each containing:
%     stim_tf (complex double, channels x frequencies x samples x trials) - Morlet wavelet
%       output of the LFP from OpenEphys_BaseAnalysis*.m (saved -v7.3, chunk layout
%       [allChannels 1 1 allTrials]). 1 kHz, 6001 samples = LFP samples 2000:8000, so the
%       stimulus (LFP sample 5001) is TF sample 3002 (stim_onset) and the last sample is
%       +2999 ms (the "times = -3000:3000" axis in BaseAnalysis is 1 ms off).
%       Frequencies are not saved in the file: rebuilt from the number of frequency rows -
%       151 rows = linspace(1,150,151) (BaseAnalysis default; steps of 0.993 Hz, so row
%       number is not exactly Hz), 150 rows = 1:150; any other count stops the script
%       (custom range in BaseAnalysis - frequencies unknown). Channel index is assumed to
%       equal the raw channel ID [inferred: holds when ProbeInfo.ChanIds = 1:Chans]
%   For tf_recompute, the same recording's LFP file, found by replacing \TF\ with \LFP\ and
%   _TF_results.mat with -LFP.mat in the TF file's path, containing:
%     stim_lfp_stimchunks (double, trials x samples x channels) - all trials of the recording
%       (same count as tr_keep); samples 2000:8000 are the TF input. The recomputation assumes
%       BaseAnalysis's default wavelet (Morlet, 4 s long, cycles logspace 6-10 over the file's
%       frequency rows, all trials joined end to end, single precision) - verified per file
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
%   stim_onset, P1wind, P2wind, Allwind, basewind (double) - stimulus sample (3002) and the
%     windows as samples: 75-350, 350-2999 and 75-2999 ms post stimulus; baseline
%     -2102..-102 ms (= samples 900:2900, as in the old TF script)
%   frex (double, 1 x nFreq) - frequency of each row of stim_tf, per file
%   bandpow (cell, 1 x nBands) - per file, mean power over the band's rows,
%     kept channels x samples x good trials
%   metriccols (struct) - one field per measure (column name), each a column over all rows
%
% Outputs:
%   tableres (table, nRows x (nMeasures + 11 + nLabels + nScores)) - one row per good
%   trial x kept channel, ordered condition -> file -> channel -> trial. Columns:
%     (measure columns, set in the calculation block; names carry "TF" so LFP/MUA/TF
%     tables can be combined:)
%     P1_TF<band>_Pow, P2_TF<band>_Pow, All_TF<band>_Pow (double, 3 per tf_bands row) -
%       mean power over the band's frequency rows and the P1 (75-350 ms post stimulus),
%       P2 (350-2999 ms) or All (75-2999 ms) window, divided by the same trial's mean band
%       power over the baseline (-2102..-102 ms). Ratio, unitless; 1 = no change
%       (e.g. P1_TFAlphaBeta_Pow, All_TFLoloGamma_Pow)
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
% Dependencies: OpenEphys_BaseAnalysis*.m (writes *_TF_results.mat, *-LFP.mat and <animal>-ProbeInfo.mat);
%   OpenEphys_editProbeInfo_ChRemove.m (sets Ch_Remove);
%   ConditionalTrialRemove.m / ConditionalTrialRemove_callable.m (tr_remove mask format,
%   tr_remove_conditional labels; TF files with an empty tr_remove must be migrated by it
%   first); ProbeInfo.ShankLabels from the channel/shank labelling scripts.
%   Requires uipickfiles (File Exchange).

%list mice to be analyzed. Note the order as it will be treated as a factor (R style)
%names must match the start of the data file names (animal-stim_TF_results.mat)
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
datadir = 'E:\Roy\Processed Silicon Probe Data\TF'; %where the result files are picked from
datavar = 'stim_tf';                                 %data variable in each file: channels x frequencies x samples x trials (complex)
tf_bands = {"AlphaBeta",[10 18];                     %band name (goes into the column names), [low high] Hz inclusive
            "LoloGamma",[30 50]};                    %add rows for more bands
tf_recompute = true;  %recompute the band wavelet output from the recording's LFP file (seconds per file, instead
                      %of minutes reading stim_tf); checked against stim_tf per file, falls back to reading it

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

stim_onset = 3002; %TF sample (1 ms) of stimulus onset: TF trials are LFP samples 2000:8000, stimulus = LFP sample 5001
P1wind = stim_onset + (75:350);        %ms post stimulus 75-350
P2wind = stim_onset + (350:2999);      %350-2999 (last TF sample is +2999 ms)
Allwind = stim_onset + (75:2999);      %75-2999
basewind = stim_onset + (-2102:-102);  %baseline, = samples 900:2900 as in the old TF script
bandsprinted = false;                  %print the rows/Hz used for each band once

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

        % ============ DATA-TYPE-SPECIFIC LOAD ============
        %frequency of each row of stim_tf (not saved in the file): rebuilt from the row count
        nfreq = size(dat,datavar,2);
        nsamp = size(dat,datavar,3);
        assert(nsamp >= Allwind(end), '%s: stim_tf has %d samples per trial, the windows need %d (6001 = LFP samples 2000:8000)', ...
            fname, nsamp, Allwind(end))
        if nfreq == 151
            frex = linspace(1,150,151); %OpenEphys_BaseAnalysis*.m default (steps of 0.993 Hz)
        elseif nfreq == 150
            frex = 1:150;
        else
            error('%s: stim_tf has %d frequency rows - expected 151 (BaseAnalysis default) or 150; frequencies unknown', fname, nfreq)
        end

        %frequency rows of each band
        bandrows = cell(1,size(tf_bands,1));
        for b = 1:size(tf_bands,1)
            bandrows{b} = find(frex >= tf_bands{b,2}(1) & frex <= tf_bands{b,2}(2));
            assert(~isempty(bandrows{b}), '%s: no frequency rows in band %s [%g %g] Hz', fname, tf_bands{b,1}, tf_bands{b,2})
            if ~bandsprinted
                fprintf('Band %s: rows %d-%d = %.2f-%.2f Hz\n', tf_bands{b,1}, bandrows{b}(1), bandrows{b}(end), ...
                    frex(bandrows{b}(1)), frex(bandrows{b}(end)))
            end
        end
        bandsprinted = true;
        %matfile can only read evenly spaced ranges: stim_tf reads use the contiguous channel/trial block, then pick
        tr_rng = min(trs):max(trs);
        ch_rng = min(chs):max(chs);

        %Fast path: recompute the band rows from the recording's LFP file exactly as OpenEphys_BaseAnalysis*.m
        %builds stim_tf (Morlet wavelets, default cycles 6-10, all trials of the file joined end to end,
        %single precision). Checked against 3 stored samples of stim_tf; if they differ (e.g. custom TF
        %settings were used), or the LFP file is missing, stim_tf is read instead.
        lfpfile = strrep(strrep(superCondis_dir{con}{file},'\TF\','\LFP\'),'_TF_results.mat','-LFP.mat');
        userecompute = tf_recompute && isfile(lfpfile);
        if tf_recompute && ~userecompute
            warning('%s: LFP file not found (%s) - reading stim_tf instead (slow)', fname, lfpfile)
        end
        if userecompute
            L = load(lfpfile,'stim_lfp_stimchunks');
            assert(size(L.stim_lfp_stimchunks,1) == numel(tr_keep), '%s: LFP file has %d trials but tr_keep has %d', ...
                fname, size(L.stim_lfp_stimchunks,1), numel(tr_keep))
            lfp = single(L.stim_lfp_stimchunks(:,2000:8000,chs)); %all trials (the wavelet sees them joined), kept channels
            clear L
            ntrall = size(lfp,1);
            wavtime = -2:1/1000:2;
            halfwave = (numel(wavtime)-1)/2;
            nConv = numel(wavtime) + nsamp*ntrall - 1;
            wavwidth = logspace(log10(6),log10(10),nfreq) ./ (2*pi*frex); %BaseAnalysis default wavelet widths
            dX = fft(reshape(permute(lfp,[2 1 3]),nsamp*ntrall,nch),nConv); %columns = channels, trials joined
            clear lfp

            %check: first row of the first band at 3 samples, all kept channels and good trials
            fi = bandrows{1}(1);
            wX = fft(exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*wavwidth(fi)^2)),nConv).';
            as = ifft((wX./max(wX)) .* dX);
            as = reshape(as(halfwave+1:end-halfwave,:),nsamp,ntrall,nch);
            chksamp = 1002:2000:5002;
            stored = dat.(datavar)(ch_rng,fi,chksamp,tr_rng);
            stored = permute(reshape(stored(chs - ch_rng(1) + 1,1,:,trs - tr_rng(1) + 1),nch,numel(chksamp),ntr),[1 3 2]);
            recomp = permute(as(chksamp,trs,:),[3 2 1]);
            reldiff = max(abs(recomp - stored),[],'all') / max(abs(stored),[],'all');
            if reldiff > 1e-4
                warning('%s: recomputed wavelet output differs from stim_tf (relative %.1e; custom TF settings?) - reading stim_tf instead (slow)', ...
                    fname, reldiff)
                userecompute = false;
            end
        end

        %mean power over each band's rows, for the kept channels and good trials: channels x samples x trials
        bandpow = cell(1,size(tf_bands,1));
        if userecompute
            fprintf('%s: band power recomputed from the LFP file (check vs stim_tf: %.1e)\n', fname, reldiff)
            for b = 1:size(tf_bands,1)
                acc = zeros(nch,nsamp,ntr);
                for fi = bandrows{b}
                    wX = fft(exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*wavwidth(fi)^2)),nConv).';
                    as = ifft((wX./max(wX)) .* dX);
                    as = reshape(as(halfwave+1:end-halfwave,:),nsamp,ntrall,nch);
                    acc = acc + double(permute(abs(as(:,trs,:)).^2,[3 1 2]));
                end
                bandpow{b} = acc / numel(bandrows{b});
            end
            clear dX as wX
        else
            %reading path: the file is stored in chunks of all channels x 1 frequency x 1 sample x all trials,
            %so all channels are read at once (same cost as one), one frequency row at a time, and turned
            %into power straight away
            fprintf('%s: reading band power from stim_tf\n', fname)
            for b = 1:size(tf_bands,1)
                acc = zeros(nch,nsamp,ntr);
                for fi = bandrows{b}
                    blk = dat.(datavar)(ch_rng,fi,:,tr_rng);
                    blk = blk(chs - ch_rng(1) + 1,1,:,trs - tr_rng(1) + 1);
                    acc = acc + reshape(abs(blk).^2,nch,nsamp,ntr);
                end
                bandpow{b} = acc / numel(bandrows{b});
            end
            clear blk
        end
        clear acc
        % =================================================

        for ch = 1:nch
            chmetrics = struct(); %one field per measure, each ntr x 1 (one value per trial)

            % ============ DATA-TYPE-SPECIFIC CALCULATION ============
            %compute the measures for this channel and set chmetrics.<ColumnName> = ntr x 1 column for each
            %(windows: P1wind / P2wind / Allwind / basewind, as samples)

            %band power in each window / the same trial's baseline band power (1 = no change)
            for b = 1:size(tf_bands,1)
                pw = reshape(permute(bandpow{b}(ch,:,:),[3 2 1]),ntr,[]); %trials x samples
                base = mean(pw(:,basewind),2);
                chmetrics.(char("P1_TF" + tf_bands{b,1} + "_Pow")) = mean(pw(:,P1wind),2) ./ base;
                chmetrics.(char("P2_TF" + tf_bands{b,1} + "_Pow")) = mean(pw(:,P2wind),2) ./ base;
                chmetrics.(char("All_TF" + tf_bands{b,1} + "_Pow")) = mean(pw(:,Allwind),2) ./ base;
            end
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
        clear bandpow pw
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
