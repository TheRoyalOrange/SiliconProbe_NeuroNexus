% NNx_StochRes_MUA_singleanimal.m
%
% Description: first look for signs of stochastic resonance (SR) in V1 multi-unit activity (MUA), one
%   animal at a time. Compares two condition groups (e.g. light alone vs light + whisker, where the extra
%   input could act as the "noise"), each made of one or more conditions (subgroups, e.g. intensities):
%     1) MUA (0/1 spike onsets, 30 kHz) -> 1 ms spike counts -> Gaussian-smoothed spike rate per trial
%     1b) response trials: a trial counts as a response when its rate rise in resp_ms (vs the pre_resp_ms
%        ending at pre_end_ms) beats resp_null_prct of the same measure at fake stim times in the pre-stim baseline
%        (mean z over the responsive channels, both groups pooled so they get the same threshold). Only
%        response trials are used from 2) on; the fraction of response trials per group is reported
%        (fewer misses with the added input would itself be an SR sign)
%     2) look at the trials: stacked (one below the other) and overlaid (how well they line up as is)
%     3) spike-rate spectrum (PSD) of each trial on the raw 1 ms counts, divided by the trial's mean rate
%        (a random/Poisson spike train is flat at 1, so rhythm = bumps above 1), condition means in the
%        probe layout. Total (evoked + induced) and induced only (condition mean trace taken out)
%     4) line the trials up by cross-correlation (one lag per trial, shared by all channels, capped at
%        +/- max_lag_ms), with 3 templates compared: leave-one-out mean (iterated), strongest trial,
%        first responsive trial. Same alignment on the pre-stim baseline as a null (alignment always
%        raises correlation a bit, even on noise)
%     5) PSD and trial-to-trial correlation (as in NNx_LFPCorrelation_LvsLW.m) on the aligned trials:
%        does one group have more consistent trials once onset jitter is taken out?
%     6) S1 (probe 2) LFP spectra next to the V1 spike-rate spectra, when you say S1 is usable (do_S1):
%        if the whisker input acts as SR "noise", V1 spiking should get more rhythmic in the bands where
%        S1 gets more power, at the same times. For group 2 (e.g. LW), sliding-window spectra of V1 spike
%        rhythmicity and of every S1 channel's power per trial; at each frequency the shapes of the two time
%        courses are correlated over a sliding block, per trial, averaged, tested against trial shuffles
%        with probe-wide cluster correction (S1-probe layout of time-frequency maps)
%        (spike-LFP / LFP-LFP phase locking was tried and dropped, see commit 1747056)
%   Sections are being built one at a time; the list above is the plan.
%
% Inputs:
%   User-set config (Setup section): animal, region, group names + condition lists, windows (ms),
%     smoothing, PSD, alignment and LFP settings, toggles
%   Spiking files picked per condition with uipickfiles,
%   <spike_dir>\<animal>\<animal>-<stim>-spiking_results.mat (saved -v7.3, read with matfile):
%     stim_spike_stimchunks (double 0/1, trials x samples x channels) - 1 at the first sample of each
%       threshold crossing, 30 kHz, 450000 samples per trial (-5 s to +10 s)
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1)
%     tr_remove (0/1, 1 x nTrials) - mask over tr_keep; only tr_keep(tr_remove == 0) are used
%     tr_remove_conditional (table, optional) - Name / TrialIdx named trial sets (NNx_TrialRemove_Conditional.m)
%   <lfp_dir>\<animal>\<animal>-<stim>-LFP.mat (do_S1 = 1, S1 channels only): stim_lfp_stimchunks (double, trials x
%     samples x channels) - peri-stimulus LFP, 1 kHz, uV, plus the same tr_keep / tr_remove
%   <probeinfo_dir>\<animal>-ProbeInfo.mat, ProbeInfo fields:
%     .Areas (cell of char, 1 x probes) - region per probe; region picks the probe
%     .ProbeMaps (cell, 1 x probes) - [rows x shanks] raw channel IDs, 0 = no site, row 1 = top
%     .Ch_Remove (cell, 1 x probes) - raw channel IDs to leave out
%     .ShankLabels (cell, 1 x probes, do_S1 = 1) - struct per probe; .(s1_shank_label) (1 x shanks) scores,
%       the S1 shank with the highest score gives the S1 reference channel (NNx_editProbeInfo_ChLabel_*.m)
%   INFERRED DATA CONTRACTS (not stated in the files, same as the other scripts that read them):
%     - stim onset is sample 150001 of stim_spike_stimchunks (30 kHz, pre = 5 s) and sample 5001 of
%       stim_lfp_stimchunks (1 kHz)
%     - channel index (dim 3) = raw channel ID (true when ChanIds = 1:Chans), both probes in one file
%     - the spiking and LFP files of a stim share tr_keep / tr_remove (checked when the LFP is loaded)
%
% Outputs (workspace; saved in the last section; more added as sections are built):
%   prb (double) - probe used; pmap (double, rows x shanks) - its ProbeMaps raw channel IDs
%   chanlist (double, channels x 1) - kept raw channel IDs (probe minus Ch_Remove), sorted
%   chrow, chshank (double, channels x 1) - ProbeMaps row (1 = top) and shank of each kept channel
%   condis (cell, 1 x conditions) - condition names, group 1's first; condgrp (1 x conditions) - group 1 or 2
%   t_bin (double, 1 x bins) - ms (vs stim onset) where each count bin starts
%   spkcount (cell, 1 x conditions) - (single, trials x bins x channels) spike counts per bin_ms bin,
%     good trials only, files of a condition pooled, channels in chanlist order
%   trinfo (cell, 1 x conditions) - (double, trials x 2) file number within the condition, trial (tr_keep value)
%   ntrials (double, 1 x conditions) - good trials per condition
%   rate (cell, 1 x conditions) - (single, trials x bins x channels) spikes/s, spkcount smoothed with a
%     Gaussian kernel smooth_ms wide (total width)
%   filesused (string, files x 4) - condition, file name, good trials, trials left out by tr_conditional_use
%   rise (cell, 1 x conditions) - (single, trials x channels) spikes/s, rate in resp_ms minus the
%     pre_resp_ms ending at pre_end_ms (from the raw counts)
%   chan_medz (single, channels x 1) - each channel's median trial z (rise vs its fake no-stim rises,
%     z = (rise - mean) / SD of the fake rises)
%   respchans (logical, channels x 1) - true = responsive channel (chan_medz >= chan_z_min)
%   respscore (cell, 1 x conditions) - (double, trials x 1) mean z over the responsive channels
%   resp (cell, 1 x conditions) - (logical, trials x 1) true = response trial (respscore > respthr).
%     Only these trials are used from here on
%   respnull (double, fake windows x 1) - the same score at the fake no-stim times (both groups pooled);
%     respthr (double) - its resp_null_prct percentile
%   respfrac (double, 1 x conditions), respfrac_grp (double, 1 x groups) - fraction of response trials
%   rate_mean, rate_sem (cell, 1 x conditions) - (single, bins x channels) spikes/s, mean and SEM of
%     rate over the response trials
%   psd_f (double, freqs x 1) - Hz of the spectra (psd_hz, multitaper nfft steps)
%   psd_tot, psd_ind, psd_pre (cell, 1 x conditions) - (single, response trials x freqs x channels)
%     multitaper spectrum of the raw counts in psd_ms (total; induced = minus the leave-one-out mean of the
%     other response trials; pre = same-length window before the stimulus), divided by the trial's
%     random (Poisson) firing level 2*mean count/fs: 1 = no rhythm. NaN at line_hz if line_out = 1
%   psdm_tot, psdm_ind, psdm_pre (cell, 1 x conditions) - (double, freqs x channels) mean over trials
%   psdm_evk (cell, 1 x conditions) - (double, freqs x channels) spectrum of the condition mean trace
%     (evoked / phase-locked part), same scale as psdm_tot (random-firing floor ~1/trials)
%   (do_align = 1) align_methods (cell, 1 x 3) - 'loo' (Woody, mean of the other trials), 'strongest'
%     (highest respscore trial), 'first' (first response trial) templates
%   lags, lags_null (cell, 1 x methods){1 x conditions} - (double, response trials x 1) ms, one shift per
%     trial from cross-correlating the responsive-channel summed rate in align_ms (null: same-length
%     pre-stim window) with the template, |lag| <= max_lag_ms. Aligned trial(t) = trial(t + lag)
%   align_use_ms (double, 1 x 2) - ms, the alignment window actually used: align_ms trimmed by max_lag_ms
%     where it would run off the loaded data (load_ms); the null window is the same length, ending at 0
%   lagpk (cell, 1 x methods){1 x conditions} - (double, response trials x 1) peak cross-correlation r
%   lagsd, lagsd_null (double, methods x conditions) - ms, SD of the lags (trial-to-trial jitter)
%   (do_align = 1, aligned-trial section) al_names (cell, 1 x versions) - 'unaligned', 'strongest'
%   psdal_f (double, freqs x 1) - Hz; psdal_evk, psdal_ind (cell, 1 x versions){1 x conditions} -
%     (double, freqs x channels) evoked / induced spectrum of the raw counts in align_use_ms, each trial
%     read at t + its lag, normalized as psdm_* (1 = random firing)
%   tc_cond, tc_cond_null (cell, 1 x versions) - (double, channels x conditions) mean pairwise trial r
%     (smoothed rate time courses, Fisher z mean, back to r) within each condition; null = pre-stim window
%   tc_trial (cell, 1 x versions){1 x conditions} - (double, response trials x channels) each trial's mean r
%     with the other trials of its condition
%   tc_grp, tc_grp_null (cell, 1 x versions) - (double, channels x groups) the same within each group
%   tc_p (cell, 1 x versions) - (double, channels x 1) label-shuffle permutation p, group 2 vs group 1
%   (do_S1 = 1) lfp_s1 (cell, 1 x conditions) - (single, response trials x samples x S1 channels) uV, S1 LFP
%     in lfp_load_ms (1 kHz), same trials/order as the MUA; t_lfp (double, 1 x samples) - ms of each sample
%   chanlist_s1 (double, channels x 1) - kept S1 raw channel IDs; prb_s1 - S1 probe number
%   s1_refch (double) - raw channel ID of the S1 reference channel (kept channel s1_ref_row of the peak
%     s1_shank_label shank)
%   (do_S1 = 1, time-frequency, group 2's response trials only)
%   tf_t (double, 1 x windows) - ms, sliding-window centres used (tf_ms, trimmed to the loaded data);
%     tf_f (double, freqs x 1) - Hz (psd_hz, steps of 1000/nfft of tf_win_ms)
%   v1tr (double, freqs x windows x trials) - V1 spike rhythmicity per group 2 response trial: multitaper
%     power of the raw counts per window / that window's random-firing level, responsive-channel mean, dB
%   ntr2 (double) - group 2 response trials used
%   r_shape_ch (double, freqs x windows x S1 channels) - shape match: per trial, r over time (a shape_block_ms
%     block of windows centred on each window) between V1 rhythmicity and that S1 channel's power at each
%     frequency, Fisher-z averaged over group 2's trials (NaN where the block runs off either end). Plotted
%   shape_zpix, shape_ppix (double, freqs x windows x S1 channels) - each pixel's z and two-sided p against
%     n_shuf trial shuffles (V1 of one trial paired with another trial's S1)
%   shape_sig (logical, freqs x windows x S1 channels) - pixels in clusters (|z| past shuf_p, positive and
%     negative apart) whose summed |z| beats shape_clustthr = the (1 - clust_p) percentile of the biggest
%     chance cluster per shuffle over all S1 channels
%   run_notes (string) - things worth knowing, for the report
%
% Dependencies: NNx_BaseAnalysis*.m (spiking + LFP files, ProbeInfo); NNx_TrialRemove_Conditional(_callable).m
%   (tr_remove mask format, optional conditional sets); uipickfiles.m; Signal Processing Toolbox (pmtm,
%   xcorr); Statistics Toolbox (corr, prctile); MATLAB R2022a+ (clim, tiledlayout 'tight', boxchart)
%

%% Setup: Specify Animal, Conditions, and Analysis Info

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Who is the data? (exactly how the file names start: <animal>-<stim>-spiking_results.mat)
animal = '20260423-p12';

%which probe has the MUA you want? (an entry of ProbeInfo.Areas, or a probe number like "1")
region = "V1";

%what are the two groups, and which conditions go in each? (you'll pick each condition's file next)
%one condition per group is fine. More than one = subgroups (e.g. intensities), kept apart in the
%results and also pooled into the group
grpnames = {'L','LW'};
condis_grp{1} = {'L_8'}; %conditions in group 1
condis_grp{2} = {'LW_8'}; %conditions in group 2

%leave out trials flagged under these conditional trial set names? (made with NNx_TrialRemove_Conditional.m,
%e.g. {"whisker_twitch_artifact"}). {} = no, only tr_remove
tr_conditional_use = {};

%where are the spiking, LFP and ProbeInfo files?
spike_dir = 'E:\Roy\Processed Silicon Probe Data\Spiking';
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%Where should I, Matlab, put the results? (Desktop while we're testing)
save_directory = 'C:\Users\rmissa\Desktop\StochRes_test'; %results in \StochRes\<animal>, figures in \Figures_StochRes\<animal>
%save_directory = 'E:\Roy\Processed Silicon Probe Data';

%%%Analysis details%%%
stim_onset_30k = 150001; %sample of stimulus onset (0 ms) in stim_spike_stimchunks, 30 kHz
stim_onset_lfp = 5001; %sample of stimulus onset (0 ms) in stim_lfp_stimchunks, 1 kHz
fs_mua = 30000; %Hz, MUA sampling rate
win_ms = [0 3000]; %ms post stimulus, the response window everything is measured in
base_ms = [-3000 0]; %ms, pre-stim window, same length as win_ms (spectrum baseline + alignment null)
bin_ms = 1; %ms, spike count bin size (1 ms -> spike rate and spectra at 1 kHz)
smooth_ms = 25; %ms, total width of the Gaussian kernel for the spike-rate traces (not the SD; ~ SD x 5)

%%%Response trial details%%%
%only trials with a response are used from the trial plots on (no-response trials just muddy things).
%a trial's score = rate in resp_ms minus rate in the pre_resp_ms ending at pre_end_ms, as a z against the same
%measure at fake stim times in the pre-stim baseline (fake_ms), averaged over the responsive channels
%(tried and dropped: resp_ms [75 500], [500 1000], [500 1500], [75 300] with the "before" stretch ending at
%70 ms; pre_resp_ms 500; median/MAD z. The first approach, below, worked best. fake_ms has to move with them, see its note)
resp_ms = [0 300]; %ms post stimulus, where a trial's response is measured
pre_resp_ms = 200; %ms, length of the "before" stretch the response is measured from
pre_end_ms = 0; %ms vs the (real or fake) stimulus, where the "before" stretch ends (0 -> -200 to 0 ms, just before the stimulus)
fake_ms = -2700:300:-300; %ms, fake stimulus times in the pre-stim baseline (no-stim reference, 9 per trial). Needs fake_ms(1) + pre_end_ms - pre_resp_ms >= base_ms(1) and fake_ms(end) + resp_ms(2) <= 0
chan_z_min = 1; %a channel counts as responsive when its typical (median) trial z is at least this
resp_null_prct = 99; %a trial is a response trial when its score beats this percentile of the fake scores

%%%Spectrum details%%%
%psd_seg_ms = 1000; %commented out: Welch dropped, only ~4 segments per trial in psd_ms. Multitaper instead
psd_ms = [350 3000]; %ms post stimulus, window for the spike-rate spectrum (sustained response, after the onset burst; = LFPCorr P2)
psd_nw = 4; %multitaper time-halfbandwidth (smooths each spectrum over ~ +/- psd_nw/window length, ~1.5 Hz for 2.65 s)
psd_hz = [3 80]; %Hz, range to look at (from 3: below that the window only holds a few cycles)
%leave out line noise from everything spectral? (only needed if 50 Hz mains shows up in the data)
line_out = 0;
    % 0 = no, keep all frequencies
    % 1 = yes, leave out line_hz
line_hz = [48 52]; %Hz, 50 Hz mains band (only used when line_out = 1). One row per band, e.g. [48 52; 98 102]

%%%Trial alignment details%%%
%line the trials up by cross-correlation and redo the spectra + trial-to-trial correlation?
do_align = 1;
    % 0 = no, unaligned trials only
    % 1 = yes, with all 3 templates (leave-one-out mean, strongest trial, first responsive trial)
align_ms = [350 3000]; %ms post stimulus, the stretch the cross-correlation looks at (sustained/rhythmic part)
max_lag_ms = 50; %ms, biggest shift allowed for a trial (either way)
woody_iter = 3; %rounds of re-aligning to the leave-one-out mean (Woody's method)
%(strongest trial = highest response score, first responsive trial = first response trial, see above)
n_perm = 1000; %label shuffles for the group 1 vs group 2 trial-to-trial correlation test

%%%S1 details%%%
%(spike-LFP / LFP-LFP phase locking was dropped, the code is in commit 1747056)
% do_lfpcoh = 1; %commented out with the phase locking, the S1 comparison only needs do_S1
% ppc_hz = 2:2:80; ppc_cycles = [3 5]; ent_nbins = 18; %commented out: phase-locking settings
%compare S1 (probe 2) LFP spectra with the V1 spike-rate spectra? only turn on for animals where probe 2 was really in S1
do_S1 = 1;
    % 0 = no
    % 1 = yes
s1_region = "S1"; %ProbeInfo.Areas name of the S1 probe
lfp_load_ms = [-5000 5000]; %ms, S1 LFP stretch kept (has to hold the sliding windows and the baseline windows)
%which S1 channel for the single-channel S1 comparison (next to the S1 average)? kept channel number s1_ref_row of
%the S1 shank with the highest score for this shank label (ProbeInfo.ShankLabels, the name you gave it when labelling)
s1_shank_label = "Whisker";
s1_ref_row = 3; %which KEPT channel on that shank, counted from the top (2 = second kept channel; removed channels skipped)
% s1_ms = [150 2000]; %commented out: fixed S1 window replaced by the sliding-window (time-frequency) comparison
%sliding-window time-frequency comparison of V1 spike rhythmicity and S1 power (group 2 only, e.g. LW)
tf_win_ms = 500; %ms, sliding window length (500 -> ~2 Hz frequency steps)
tf_step_ms = 50; %ms, step between windows
tf_ms = [-1000 3000]; %ms, first and last window centre (trimmed to what the loaded MUA / LFP can hold)
tf_nw = 2; %multitaper time-halfbandwidth for the sliding windows (~ +/- tf_nw/window length = +/- 4 Hz for 500 ms)
% tf_base_ms = [-2500 -250]; %commented out: no S1 baseline needed, the shape correlation ignores offsets
shape_block_ms = 1000; %ms of window centres per shape correlation (V1 vs S1 time course at each frequency, per trial, sliding)
n_shuf = 200; %trial shuffles for the shape-map null (V1 of one trial paired with S1 of another; 200 ~ 3-6 min for ~64 S1 channels)
shuf_p = 0.05; %pixel threshold (two-sided) for building clusters
clust_p = 0.05; %cluster threshold: beat the (1 - clust_p) percentile of the biggest chance cluster

%%%Plotting details%%%
%stacked + overlaid single trials?
plot_trials = 1;
    % 0 = no
    % 1 = yes (probe layout, one figure per condition)
stack_offset = 100; %spikes/s, how far each trial is moved down from the one before in the stacked plot
plot_ms = [-500 2000]; %ms, time range shown in the trial plots
%probe-layout spectra?
plot_psd = 1;
%alignment method comparison?
plot_align = 1;
%trial-to-trial correlation maps?
plot_corr = 1;
%S1 LFP spectra vs V1 spike-rate spectra? (do_S1 = 1)
plot_lfpcoh = 1;
%save the figures?
plotsave = 0;

%delete big in-between variables as we go? (saves memory)
mr_clean = 1;

%%%nothing else to fill in%%%




%Now beginning the data loading section

%% Data Extraction: Load ProbeInfo and pick the probe + channels

run_notes = strings(0,1); %things worth knowing, shown in the report at the end (nothing printed mid-script)

%load the ProbeInfo
pf = [probeinfo_dir '\' animal '-ProbeInfo.mat'];
if isfile(pf) == false
    error(['No ProbeInfo file for ' animal ': ' pf])
end
load(pf,'ProbeInfo');

%which probe? (region name, or a probe number)
prb = find(strcmp(string(ProbeInfo.Areas), region));
if isempty(prb) %region can also be a probe number
    prb = str2double(region);
end
if numel(prb) ~= 1 || isnan(prb)
    error(['region "' char(region) '" does not pick exactly one probe in ProbeInfo.Areas'])
end

%kept channels (probe minus Ch_Remove) and where they sit on the probe
pmap = ProbeInfo.ProbeMaps{prb}; % [rows x shanks] raw channel IDs, 0 = no site, row 1 = top
chremove = [];
if numel(ProbeInfo.Ch_Remove) >= prb
    chremove = ProbeInfo.Ch_Remove{prb};
end
chanlist = setdiff(nonzeros(pmap), chremove(:)); %sorted column of kept raw channel IDs
nch = numel(chanlist);
chrow = zeros(nch,1); %ProbeMaps row of each kept channel (1 = top)
chshank = zeros(nch,1); %ProbeMaps column (shank)
for ch = 1:nch
    [chrow(ch), chshank(ch)] = find(pmap == chanlist(ch)); %where this channel sits on the probe
end


%% Data Extraction: Pick the spiking files of each condition

%line the conditions up in one list, and remember which group each is in
condis = [condis_grp{:}]; %{1 x conditions} names, group 1's first
condgrp = []; %group of each condition (1 or 2)
for g = 1:numel(condis_grp)
    condgrp = [condgrp g*ones(1,numel(condis_grp{g}))];
end
ncond = numel(condis);

%pick the file(s) for each condition (more than one file = trials pooled)
condi_files = cell(1,ncond);
for con = 1:ncond
    condi_files{con} = uipickfiles('FilterSpec',[spike_dir '\' animal],'Prompt',['Choose ' condis{con} ' spiking file(s)'],'REFilter','-spiking_results\.mat$');
    if iscell(condi_files{con}) == false
        error(['No files picked for ' condis{con}])
    end
    for file = 1:numel(condi_files{con})
        [~,fname] = fileparts(condi_files{con}{file});
        if startsWith(fname, [animal '-']) == false
            error(['Condition ' condis{con} ': file "' fname '" is not from ' animal])
        end
    end
end


%% Data Extraction: Load the good trials of each file as 1 ms spike counts

%stretch to load: the pre-stim baseline (spectra, response trials, alignment null) to the end of the
%response window, and whatever the trial plots show
load_ms = [min([base_ms(1) plot_ms(1)]) max([win_ms(2) plot_ms(2)])]; %ms around stim onset
%load_ms = [min([base_ms(1) plot_ms(1)])-max_lag_ms, max([win_ms(2) plot_ms(2)])+max_lag_ms]; %commented out: padding meant reloading for every max_lag_ms change. The alignment trims its window to fit instead
perms = fs_mua/1000; %MUA samples per ms (30)
perbin = bin_ms*perms; %MUA samples per count bin
nbins = diff(load_ms)/bin_ms; %count bins in the loaded stretch
loadsamp = (stim_onset_30k + load_ms(1)*perms) : (stim_onset_30k + load_ms(2)*perms - 1); %MUA samples pulled from each file
t_bin = load_ms(1) + (0:nbins-1)*bin_ms; %ms of the start of each bin (bin at 0 ms starts at stim onset)

spkcount = cell(1,ncond); %{conditions} [trials x bins x kept channels] spike counts per bin, single
trinfo = cell(1,ncond); %{conditions} [trials x 2] file number (within the condition), trial number (tr_keep value)
filesused = strings(0,4); %condition, file, good trials, trials left out by tr_conditional_use (for the report)

tic
for con = 1:ncond
    for file = 1:numel(condi_files{con})
        [~,fname] = fileparts(condi_files{con}{file});
        disp(['MUA: Loading ' condis{con} ', ' fname])

        sdat = matfile(condi_files{con}{file}); %-v7.3 file, read only the parts we need
        tr_keep = sdat.tr_keep;
        tr_remove = sdat.tr_remove;
        if numel(tr_remove) ~= numel(tr_keep)
            error([fname ': tr_remove is not a 0/1 mask over tr_keep (old format). Run the file through NNx_TrialRemove_Conditional first'])
        end
        szspk = size(sdat,'stim_spike_stimchunks'); % [trials samples channels] without loading it
        nsamp_file = szspk(2);
        if loadsamp(1) < 1 || loadsamp(end) > nsamp_file
            error([fname ': the windows run off the end of the trials (' num2str(nsamp_file) ' samples)'])
        end

        trs = tr_keep(tr_remove == 0); %good trials
        trs = trs(:);

        %leave out the trials flagged under tr_conditional_use (TrialIdx = tr_keep values, NaN = evaluated, none flagged)
        nexcl = 0;
        if isempty(tr_conditional_use) == false
            if any(strcmp(who('-file',condi_files{con}{file}),'tr_remove_conditional'))
                tc = sdat.tr_remove_conditional;
            else
                tc = table('Size',[0,2],'VariableTypes',{'string','double'},'VariableNames',{'Name','TrialIdx'});
            end
            for k = 1:numel(tr_conditional_use)
                namerows = strcmp(tc.Name, tr_conditional_use{k});
                if any(namerows) == false
                    run_notes(end+1) = string(fname) + ": never evaluated for " + tr_conditional_use{k} + ", none of its trials left out";
                    continue
                end
                flagged = ismember(trs, tc.TrialIdx(namerows));
                nexcl = nexcl + sum(flagged);
                trs = trs(flagged == false);
            end
        end
        ntr = numel(trs);

        %matfile only reads evenly spaced ranges, so read the block of trials + channels spanning the good
        %ones, then keep the good ones in memory.
        %the file is stored (compressed) in little pieces of ALL trials x ALL channels x 2 samples, and a
        %piece has to be unpacked whole every time any part of it is read. So read every channel at once,
        %a short stretch of time at a time (all of it at once would be GBs), instead of channel by channel
        %(that unpacked the whole file once per channel, ~60x slower)
        % for ch = 1:nch %commented out: one channel at a time re-unpacks every piece for each channel, very slow
        %     chspikes = sdat.stim_spike_stimchunks(tr_rng,loadsamp,chanlist(ch)); % [trial block x samples] 0/1 spike onsets
        %     chspikes = single(chspikes(trs - tr_rng(1) + 1,:)); %just the good trials
        %     counts(:,:,ch) = reshape(sum(reshape(chspikes,ntr,perbin,nbins),2),ntr,nbins); %spikes per bin: group each bin's samples, add them up
        % end
        tr_rng = min(trs):max(trs);
        ch_rng = min(chanlist):max(chanlist);
        blk_bins = 200/bin_ms; %bins per read (200 ms -> ~150-400 MB per read as double)
        counts = zeros(ntr,nbins,nch,'single'); % [good trials x bins x kept channels], single halves the memory, counts are small whole numbers
        for b1 = 1:blk_bins:nbins
            bidx = b1:min(b1+blk_bins-1, nbins); %bins in this read
            sidx = loadsamp((bidx(1)-1)*perbin+1 : bidx(end)*perbin); %their MUA samples
            blk = sdat.stim_spike_stimchunks(tr_rng,sidx,ch_rng); % [trial block x samples x channel block] 0/1 spike onsets
            blk = single(blk(trs - tr_rng(1) + 1, :, chanlist - ch_rng(1) + 1)); %just the good trials and kept channels
            counts(:,bidx,:) = reshape(sum(reshape(blk,ntr,perbin,numel(bidx),nch),2),ntr,numel(bidx),nch); %spikes per bin: group each bin's samples, add them up
        end
        clear blk

        %pool the files of the same condition
        spkcount{con} = cat(1, spkcount{con}, counts);
        trinfo{con} = cat(1, trinfo{con}, [file*ones(ntr,1) trs]);
        filesused(end+1,:) = [string(condis{con}) string(fname) string(ntr) string(nexcl)];
    end
end
clear counts sdat
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);

ntrials = zeros(1,ncond); %good trials per condition
for con = 1:ncond
    ntrials(con) = size(spkcount{con},1);
end
disp('MUA: All loaded, now you have 1 ms spike counts for every good trial')




%Now beginning the spike rate section

%% Data PreProcessing: Spike counts to spike rate

disp('MUA: Smoothing spike counts into spike rates')

%Gaussian kernel over the spike COUNTS (a standard kernel density firing rate). Every spike counts,
%and the kernel is smooth so it doesn't add sharp edges (which would show up in the spectra later)
%(the NNx_BaseAnalysis way, 1 if any spike in the ms + 50 ms flat moving mean, was compared and dropped:
%it undercounts bursts, when 2+ spikes land in the same ms)
rate = cell(1,ncond); %{conditions} [trials x bins x channels] spikes/s, single
for con = 1:ncond
    rate{con} = smoothdata(spkcount{con}, 2, 'gaussian', smooth_ms/bin_ms) * (1000/bin_ms); %smooths along time for every trial + channel at once
end
%(the smoothing shrinks its window at the very ends of the loaded stretch, load_ms, which is outside
%every window we use, so no padding needed)
%the raw counts (spkcount) stay for the spectra and the response trials




%Now beginning the response trial section

%% Data Analysis: Response trials (yes/no)

%Did the trial respond at all? For each trial and channel: the spike rate in resp_ms minus the rate in
%the pre_resp_ms ending at pre_end_ms (the rise). The same rise is measured at fake stimulus times
%in the pre-stim baseline (fake_ms), where nothing happened. Each channel's rises are put in units of
%its no-stim spread (z: minus the mean, divided by the SD). Responsive channels = typical (median) trial z >= chan_z_min. Trial score =
%mean z over them. Response trial = score beats resp_null_prct of the fake windows' scores.
%Both groups are pooled for the z, the channels and the threshold, so they're judged the same way.
%Uses the raw counts (spkcount), so the smoothing kernel can't smear the response into the "before" stretch
disp('MUA: Finding the trials with a response')
tic

%do the fake windows fit? (inside the loaded stretch, and ending before the real stimulus)
if fake_ms(1) + pre_end_ms - pre_resp_ms < load_ms(1)
    error(['The first fake window starts before the loaded data (' num2str(load_ms(1)) ' ms). Move fake_ms later or base_ms earlier'])
end
if fake_ms(end) + resp_ms(2) > 0
    error('The last fake response window runs past the real stimulus. Move fake_ms earlier')
end

respidx = find(t_bin >= resp_ms(1) & t_bin < resp_ms(2)); %bins of the response window
%preidx = find(t_bin >= -pre_resp_ms & t_bin < 0); %commented out: the "before" stretch now ends at pre_end_ms, not at the stimulus
preidx = find(t_bin >= pre_end_ms - pre_resp_ms & t_bin < pre_end_ms); %bins of the "before" stretch
nfake = numel(fake_ms); %fake windows per trial
tohz = 1000/bin_ms; %mean count per bin -> spikes/s

rise = cell(1,ncond); %{conditions} [trials x channels] spikes/s, rate rise after the stimulus
risefake = cell(1,ncond); %{conditions} [(fake windows*trials) x channels] the same at the fake times
for con = 1:ncond
    ntr = ntrials(con);
    rise{con} = reshape(mean(spkcount{con}(:,respidx,:),2) - mean(spkcount{con}(:,preidx,:),2), ntr, nch) * tohz;
    risefake{con} = zeros(nfake*ntr, nch);
    for f = 1:nfake
        fresp = find(t_bin >= fake_ms(f)+resp_ms(1) & t_bin < fake_ms(f)+resp_ms(2)); %fake response window
        %fpre = find(t_bin >= fake_ms(f)-pre_resp_ms & t_bin < fake_ms(f)); %commented out: the "before" stretch now ends at pre_end_ms
        fpre = find(t_bin >= fake_ms(f)+pre_end_ms-pre_resp_ms & t_bin < fake_ms(f)+pre_end_ms); %and its "before" stretch
        rows = (f-1)*ntr + (1:ntr); %where this fake window's trials go
        risefake{con}(rows,:) = reshape(mean(spkcount{con}(:,fresp,:),2) - mean(spkcount{con}(:,fpre,:),2), ntr, nch) * tohz;
    end
end

%every condition pooled (both groups)
allrise = cat(1, rise{:}); % [all trials x channels]
allfake = cat(1, risefake{:}); % [all fake windows x channels]
fakemu = mean(allfake,1); % [1 x channels] typical rise when nothing happened
fakesd = std(allfake,0,1); % [1 x channels] its spread
% fakemu = median(allfake,1); %commented out: tried robust median/MAD z, the calls looked worse than mean/SD
% fakesd = mad(allfake,1,1) * 1.4826; %(MAD x1.4826 matches an SD for normal data)
fakesd(fakesd == 0) = NaN; %a (nearly) silent channel has no spread, can't be put in z, it drops out
zch = (allrise - fakemu) ./ fakesd; % [all trials x channels] rise in units of the no-stim spread
zfake = (allfake - fakemu) ./ fakesd; % [all fake windows x channels]

%responsive channels
chan_medz = median(zch,1,'omitnan')'; % [channels x 1] typical trial z of each channel
respchans = chan_medz >= chan_z_min; % [channels x 1] true = responsive channel
if any(respchans) == false
    error(['No channel has a typical trial z >= chan_z_min (' num2str(chan_z_min) '), so no responses can be called. Try a lower chan_z_min'])
end

%trial scores and the threshold
sc = mean(zch(:,respchans),2); % [all trials x 1]
respnull = mean(zfake(:,respchans),2); % [all fake windows x 1] the same score when nothing happened
respthr = prctile(respnull, resp_null_prct); %a response has to beat this

%back to the conditions (same order they were stacked in)
respscore = cell(1,ncond); %{conditions} [trials x 1] trial score (mean z over the responsive channels)
resp = cell(1,ncond); %{conditions} [trials x 1] true = response trial
respfrac = zeros(1,ncond); %fraction of response trials per condition
row = 0;
for con = 1:ncond
    idx = row + (1:ntrials(con));
    respscore{con} = sc(idx);
    resp{con} = sc(idx) > respthr;
    respfrac(con) = mean(resp{con});
    row = row + ntrials(con);
end

%fraction of response trials per group (its conditions pooled)
respfrac_grp = zeros(1,numel(condis_grp));
for g = 1:numel(condis_grp)
    respfrac_grp(g) = sum(cat(1, resp{condgrp == g})) / sum(ntrials(condgrp == g));
end

%condition mean +/- SEM of the spike rate, response trials only
rate_mean = cell(1,ncond); %{conditions} [bins x channels] spikes/s
rate_sem = cell(1,ncond); %{conditions} [bins x channels] spikes/s
for con = 1:ncond
    rtr = resp{con} == true; %response trials
    rate_mean{con} = reshape(mean(rate{con}(rtr,:,:),1), nbins, nch);
    rate_sem{con} = reshape(std(rate{con}(rtr,:,:),0,1), nbins, nch) / sqrt(sum(rtr));
    if sum(rtr) < 2
        run_notes(end+1) = string(condis{con}) + ": only " + sum(rtr) + " response trial(s), its mean / SEM are not worth much";
    end
end

clear allrise allfake zfake sc
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data Visualization: Response trial check
%(temporary section, to judge the response trial method. Not part of the final script)

disp('MUA: Making some pictures of which trials got called a response')

if plotsave == 1
    %create figure folder
    if isfolder([save_directory '\Figures_StochRes\' animal]) == false
        mkdir([save_directory '\Figures_StochRes\' animal])
    end
end

pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
t_plot = t_bin(pidx); %ms
grpcol = [0 0.45 0.74; 0.85 0.33 0.10]; %colour of group 1, group 2

%picture 1: trial scores vs the fake (no stim) scores
allsc = [respnull; cat(1, respscore{:})];
edges = linspace(min(allsc), max(allsc), 60); %same bins for every histogram
figure()
histogram(respnull,edges,'Normalization','probability','FaceColor',[.6 .6 .6],'EdgeColor','none','DisplayName','fake windows (no stim)')
hold on
for g = 1:numel(condis_grp)
    gsc = cat(1, respscore{condgrp == g}); %all of this group's trial scores
    histogram(gsc,edges,'Normalization','probability','FaceColor',grpcol(g,:),'FaceAlpha',.5,'EdgeColor','none', ...
        'DisplayName',[grpnames{g} ' trials (' num2str(round(100*respfrac_grp(g))) '% response)'])
end
xline(respthr,'k--','LineWidth',1.5,'DisplayName',['threshold (' num2str(resp_null_prct) 'th percentile of fake)'])
hold off
xlabel('trial score (mean z over the responsive channels)')
ylabel('fraction of trials / windows')
legend('Location','northeast')
title([animal ' ' ProbeInfo.Areas{prb} '-Probe: response trial scores vs no-stim (' num2str(sum(respchans)) ' responsive channels)'])

if plotsave == 1
    fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'resptrials_scores']);
    savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
end

%picture 2: every trial of a condition, which got called what
for con = 1:ncond
    trsum = sum(rate{con}(:,pidx,respchans),3); % [trials x bins shown] rate summed over the responsive channels
    rtr = resp{con} == true;

    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(3,12,'TileSpacing','compact','Padding','compact');

    %all trials as an image, in recording order
    ax1 = nexttile(tl,1,[2 11]);
    imagesc(ax1, t_plot, 1:ntrials(con), trsum)
    xline(ax1,0,'w')
    xline(ax1,resp_ms,'w--')
    cb = colorbar(ax1);
    cb.Label.String = 'spikes/s, summed over responsive channels';
    ylabel(ax1,'trial (recording order)')
    set(ax1,'XTickLabel',[])
    title(ax1,[animal ' ' condis{con} ': ' num2str(sum(rtr)) ' of ' num2str(ntrials(con)) ' trials called a response'])

    %strip beside it: green = response, grey = no response
    ax2 = nexttile(tl,12,[2 1]);
    imagesc(ax2, 1, 1:ntrials(con), double(rtr))
    colormap(ax2,[0.8 0.8 0.8; 0 0.6 0])
    clim(ax2,[0 1])
    set(ax2,'XTick',[],'YTick',[])
    title(ax2,'resp?')

    %mean trace of the response vs no-response trials
    ax3 = nexttile(tl,25,[1 11]);
    plot(ax3, t_plot, mean(trsum(rtr,:),1), 'Color',[0 0.6 0], 'LineWidth',1.5)
    hold(ax3,'on')
    plot(ax3, t_plot, mean(trsum(rtr == false,:),1), 'Color',[.5 .5 .5], 'LineWidth',1.5)
    xline(ax3,0)
    xline(ax3,resp_ms,'--')
    hold(ax3,'off')
    xlim(ax3,plot_ms)
    xlabel(ax3,'ms from stim onset')
    ylabel(ax3,'spikes/s')
    legend(ax3,{['response (' num2str(sum(rtr)) ')'],['no response (' num2str(sum(rtr == false)) ')']},'Location','northeast')

    if plotsave == 1
        fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-' 'resptrials_trials']);
        savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
    end
end

%picture 3: each channel's typical z on the probe, responsive channels circled
zmap = nan(size(pmap)); % [rows x shanks]
for ch = 1:nch
    zmap(chrow(ch),chshank(ch)) = chan_medz(ch);
end
figure()
h = imagesc(zmap);
set(h,'AlphaData',isnan(zmap) == false) %no site / removed = blank
cb = colorbar;
cb.Label.String = 'median trial z';
hold on
for ch = 1:nch
    text(chshank(ch), chrow(ch), num2str(chanlist(ch)), 'HorizontalAlignment','center','FontSize',8)
end
plot(chshank(respchans), chrow(respchans), 'ko', 'MarkerSize',22, 'LineWidth',1.5) %responsive channels
[rrow,rshank] = find(pmap > 0 & ismember(pmap, chanlist) == false); %removed channels
for i = 1:numel(rrow)
    text(rshank(i), rrow(i), 'x', 'HorizontalAlignment','center','Color',[.5 .5 .5])
end
hold off
set(gca,'XTick',1:size(pmap,2),'YTick',1:size(pmap,1))
xlabel('shank')
ylabel('row (1 = top)')
title([animal ' ' ProbeInfo.Areas{prb} '-Probe: typical response z per channel (o = responsive, z >= ' num2str(chan_z_min) '; x = removed)'])

if plotsave == 1
    fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'resptrials_channels']);
    savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
end



%% Data Visualization: Response trial check, every trial on one channel
%(temporary section, to judge the response trial method. Not part of the final script)

check_ch = 22; %raw channel ID to show
raw_offset = 1000/bin_ms; %spikes/s, how far each trial is moved down in the unsmoothed stacked plot (1 spike in a bin = 1000/bin_ms spikes/s)

disp(['MUA: Plotting every trial on CH' num2str(check_ch) ', by response call'])

if plotsave == 1
    %create figure folder
    if isfolder([save_directory '\Figures_StochRes\' animal]) == false
        mkdir([save_directory '\Figures_StochRes\' animal])
    end
end

chi = find(chanlist == check_ch); %its place in the channel list
if isempty(chi)
    error(['check_ch ' num2str(check_ch) ' is not one of the kept channels of this probe'])
end
pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
t_plot = t_bin(pidx); %ms

for con = 1:ncond
    rtr = resp{con} == true;

    %the two versions to show: smoothed rate, and the raw counts as a rate (no smoothing)
    versions = {rate{con}(:,pidx,chi), spkcount{con}(:,pidx,chi) * (1000/bin_ms)}; % each [trials x bins shown] spikes/s
    vnames = {['Gaussian-smoothed (' num2str(smooth_ms) ' ms)'], ['not smoothed (' num2str(bin_ms) ' ms bins)']};
    voffset = [stack_offset raw_offset]; %stacking distance for each version
    vfile = {'resptrials_lines','resptrials_lines_raw'};

    for v = 1:2
        trdat = versions{v}; % [trials x bins shown]
        off = voffset(v);

        figure()
        set(gcf, 'Position', get(0, 'Screensize'));

        %left: trials stacked, one below the other (trial 1 on top), green = response, grey = no response
        subplot(1,2,1)
        hold on
        for tr = 1:ntrials(con)
            if rtr(tr) == true
                col = [0 0.6 0]; %response
            else
                col = [.6 .6 .6]; %no response
            end
            plot(t_plot, trdat(tr,:) - (tr-1)*off, 'Color',col, 'LineWidth',1) %each trial moved down by off
        end
        xline(0)
        xline(resp_ms,'--')
        hold off
        xlim(plot_ms)
        ylim([-(ntrials(con))*off - off, max(trdat(1,:)) + off])
        set(gca,'YTick',fliplr(-(0:ntrials(con)-1)*off),'YTickLabel',fliplr(trinfo{con}(:,2)')) %label each line with its trial number
        xlabel('ms from stim onset')
        ylabel(['trial (each moved down ' num2str(off) ' spikes/s)'])
        title('stacked')

        %right: all trials on top of each other (no-response first, so the response trials sit on top)
        subplot(1,2,2)
        hold on
        for tr = find(rtr == false)'
            plot(t_plot, trdat(tr,:), 'Color',[.6 .6 .6], 'LineWidth',0.75)
        end
        for tr = find(rtr == true)'
            plot(t_plot, trdat(tr,:), 'Color',[0 0.6 0], 'LineWidth',0.75)
        end
        xline(0)
        xline(resp_ms,'--')
        hold off
        xlim(plot_ms)
        xlabel('ms from stim onset')
        ylabel('spikes/s')
        title(['overlaid: green = response (' num2str(sum(rtr)) '), grey = no response (' num2str(sum(rtr == false)) ')'])

        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' condis{con} ' CH' num2str(check_ch) ': every trial, ' vnames{v} ', by response call'])

        if plotsave == 1
            fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-CH' num2str(check_ch) '-' vfile{v}]);
            savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
        end
    end
end




%Now beginning the single trial pictures (response trials only from here on)

%% Data Visualization: Stacked trials, one channel at a time

%every response trial's spike rate, one below the other (first trial on top), one channel at a time,
%conditions side by side. Channels go from the strongest response (highest median z) down.
%Flip through with the left/right arrow keys or the slider at the bottom, s = save the channel you're on,
%Esc or close the window = done. (the script waits here until you're done)
if plot_trials == 1
    disp('MUA: Stacked trial pictures, flip through the channels (Esc or close the window when done)')

    pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
    t_plot = t_bin(pidx); %ms
    [~,chorder] = sort(chan_medz,'descend','MissingPlacement','last'); %channels, strongest response first

    %the figure: a key press or a slider move gets stored in the figure's UserData, the loop below reacts to it
    stackfig = figure();
    set(stackfig, 'Position', get(0, 'Screensize'));
    set(stackfig, 'WindowKeyPressFcn', @(src,evt) set(src,'UserData',evt.Key)); %key press -> its name (works even when the slider has focus)
    sld = uicontrol(stackfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.005 0.4 0.025], ...
        'Min',1, 'Max',max(nch,2), 'Value',1, 'SliderStep',[1 10]/max(nch-1,1), ...
        'Callback', @(src,evt) set(stackfig,'UserData',round(src.Value))); %slider move -> the rank to jump to

    k = 1; %rank of the channel shown (1 = strongest response)
    while true
        ch = chorder(k); %place in chanlist of the channel to show

        %draw it (clear the last channel's panels, the slider stays)
        set(0,'CurrentFigure',stackfig)
        delete(findobj(stackfig,'Type','axes'))
        for con = 1:ncond
            ax = subplot(1,ncond,con);
            rtrs = find(resp{con} == true); %response trials
            nresp = numel(rtrs);
            if nresp == 0
                axis(ax,'off')
                text(ax,0.5,0.5,[condis{con} ': no response trials'],'Units','normalized','HorizontalAlignment','center','Interpreter','none')
                continue
            end
            trdat = reshape(rate{con}(rtrs,pidx,ch), nresp, numel(pidx)); % [response trials x bins shown] spikes/s
            stackshift = (0:nresp-1)' * stack_offset; % [response trials x 1] how far down each trial goes
            plot(ax, t_plot, (trdat - stackshift)', 'k', 'LineWidth',0.75) %all trials in one go, each column is a trial
            hold(ax,'on')
            xline(ax,0,'r')
            hold(ax,'off')
            xlim(ax,plot_ms)
            ylim(ax,[-nresp*stack_offset, max(trdat(1,:)) + stack_offset])
            set(ax,'YTick',fliplr(-(0:nresp-1)*stack_offset),'YTickLabel',fliplr(trinfo{con}(rtrs,2)')) %each line labelled with its trial number
            xlabel(ax,'ms from stim onset')
            if con == 1
                ylabel(ax,['trial (each moved down ' num2str(stack_offset) ' spikes/s)'])
            end
            title(ax,[condis{con} ' (' num2str(nresp) ' of ' num2str(ntrials(con)) ' trials)'],'Interpreter','none')
        end
        if respchans(ch) == true
            rtxt = 'responsive';
        else
            rtxt = 'not responsive';
        end
        sgtitle(stackfig, {[animal ' ' ProbeInfo.Areas{prb} '-Probe CH' num2str(chanlist(ch)) ': rank ' num2str(k) ' of ' num2str(nch) ...
            ', median z = ' num2str(chan_medz(ch),'%.2f') ' (' rtxt ')'], ...
            'left/right arrow keys or slider = change channel,  s = save this one,  Esc or close = done'}, 'Interpreter','none')
        sld.Value = k; %keep the slider in step with the arrow keys

        %wait for the next key press / slider move / window close
        set(stackfig,'UserData','') %clear the last one, so the same key twice still counts
        waitfor(stackfig,'UserData')
        if isvalid(stackfig) == false %window closed
            break
        end
        key = get(stackfig,'UserData');
        if isnumeric(key) %slider
            k = key;
        elseif strcmp(key,'rightarrow') || strcmp(key,'downarrow')
            k = min(k+1, nch);
        elseif strcmp(key,'leftarrow') || strcmp(key,'uparrow')
            k = max(k-1, 1);
        elseif strcmp(key,'s') %save the channel you're on
            if isfolder([save_directory '\Figures_StochRes\' animal]) == false
                mkdir([save_directory '\Figures_StochRes\' animal])
            end
            fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-CH' num2str(chanlist(ch)) '-' 'trials_stacked']);
            savefig(stackfig,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
            disp(['MUA: saved ' fname])
        elseif strcmp(key,'escape')
            close(stackfig)
            break
        end
    end
    disp('MUA: Done with the stacked trials')
end


%% Data Visualization: Overlaid trials

%the response trials all on top of each other with nothing shifted (grey), plus the condition mean
%(black), every channel in the probe layout. How well do the trials line up as they are?
if plot_trials == 1
    disp('MUA: Making overlaid trial pictures')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
    t_plot = t_bin(pidx); %ms
    [nrow,nshank] = size(pmap);

    for con = 1:ncond
        rtrs = find(resp{con} == true); %response trials
        nresp = numel(rtrs);
        if nresp == 0
            run_notes(end+1) = string(condis{con}) + ": no response trials, no overlaid trial picture";
            continue
        end
        trdat = rate{con}(rtrs,pidx,:); % [response trials x bins shown x channels] spikes/s
        ytop = max(trdat,[],'all') * 1.05; %one y range for the whole probe

        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout(nrow,nshank,'TileSpacing','none','Padding','tight');
        for row = 1:nrow
            for shank = 1:nshank
                chid = pmap(row,shank); %raw channel ID at this spot
                ax = nexttile(tl,(row-1)*nshank + shank);
                if chid == 0 %no site here
                    axis(ax,'off')
                    continue
                end
                ch = find(chanlist == chid);
                if isempty(ch) %removed channel
                    axis(ax,'off')
                    text(ax,0.5,0.5,['CH' num2str(chid) ' removed'],'Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5])
                    continue
                end
                plot(ax, t_plot, reshape(trdat(:,:,ch),nresp,[])', 'Color',[.65 .65 .65], 'LineWidth',0.5) %every response trial in one go (columns = trials)
                hold(ax,'on')
                plot(ax, t_plot, rate_mean{con}(pidx,ch), 'k', 'LineWidth',1.5) %condition mean (response trials)
                xline(ax,0,'r')
                hold(ax,'off')
                xlim(ax,plot_ms)
                ylim(ax,[0 ytop])
                text(ax,0.02,0.92,['CH' num2str(chid)],'Units','normalized','Color','b')
                if row < nrow
                    set(ax,'XTickLabel',[])
                end
                if shank > 1
                    set(ax,'YTickLabel',[])
                end
            end
        end
        xlabel(tl,'ms from stim onset')
        ylabel(tl,'spikes/s')
        title(tl,[animal ' ' ProbeInfo.Areas{prb} '-Probe ' condis{con} ': overlaid response trials (grey, ' num2str(nresp) ' of ' num2str(ntrials(con)) ') and their mean (black)'],'Interpreter','none')

        if plotsave == 1
            fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-' 'trials_overlaid']);
            savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
        end
    end
end




%Now beginning the spectrum section

%% Data Analysis: Spike-rate spectrum (PSD) per trial and condition mean

%Is the firing rhythmic, and at what frequency? Multitaper spectrum (pmtm) of each response trial's raw
%1 ms spike counts in psd_ms (the sustained response, after the onset burst), per channel.
%Each trial's spectrum is divided by what random (Poisson) firing at that trial's own rate would give,
%so 1 = no rhythm, > 1 = rhythm at that frequency, and conditions with different rates compare fairly.
%Three versions per trial:
%  total   - the trial as it is
%  induced - the trial minus the mean of the OTHER response trials. Rhythm at the same phase on every
%            trial is in the mean and gets taken out; rhythm whose phase differs between trials
%            (e.g. jittered) mostly cancels in the mean, so it stays here
%  pre     - a same-length stretch just before the stimulus (no-stim reference)
%plus the evoked spectrum: the spectrum of the condition mean trace itself (the phase-locked part),
%on the same scale as total (its random-firing floor is ~1/trials instead of 1, since averaging cancels noise)
disp('MUA: Spike-rate spectra')
tic

fs = 1000/bin_ms; %Hz, sampling rate of the counts
psdidx = find(t_bin >= psd_ms(1) & t_bin < psd_ms(2)); %bins of the response window
psdpre_ms = [-diff(psd_ms) 0]; %ms, same-length window just before the stimulus
if psdpre_ms(1) < load_ms(1)
    error(['The pre-stim spectrum window (' num2str(psdpre_ms(1)) ' ms) starts before the loaded data (' num2str(load_ms(1)) ' ms). Shorten psd_ms or load more with base_ms'])
end
psdpreidx = find(t_bin >= psdpre_ms(1) & t_bin < psdpre_ms(2)); %bins of the pre-stim window
nwin = numel(psdidx); %samples per window
nfft = 2^nextpow2(nwin);

%frequencies kept (pmtm gives 0 to fs/2 in nfft/2+1 steps)
f_all = (0:nfft/2)' * fs/nfft; %Hz
fkeep = f_all >= psd_hz(1) & f_all <= psd_hz(2);
psd_f = f_all(fkeep); % [freqs x 1] Hz
nf = numel(psd_f);
linemask = false(nf,1); %line noise bins, set to NaN if line_out = 1
if line_out == 1
    for b = 1:size(line_hz,1)
        linemask(psd_f >= line_hz(b,1) & psd_f <= line_hz(b,2)) = true;
    end
end

psd_tot = cell(1,ncond); %{conditions} [response trials x freqs x channels] total, normalized (1 = random firing)
psd_ind = cell(1,ncond); %{conditions} [response trials x freqs x channels] induced, normalized
psd_pre = cell(1,ncond); %{conditions} [response trials x freqs x channels] pre-stim, normalized
psdm_tot = cell(1,ncond); %{conditions} [freqs x channels] mean over trials
psdm_ind = cell(1,ncond);
psdm_pre = cell(1,ncond);
psdm_evk = cell(1,ncond); %{conditions} [freqs x channels] spectrum of the condition mean trace (evoked)
for con = 1:ncond
    rtrs = find(resp{con} == true); %response trials
    nresp = numel(rtrs);
    psd_tot{con} = nan(nresp,nf,nch,'single');
    psd_ind{con} = nan(nresp,nf,nch,'single');
    psd_pre{con} = nan(nresp,nf,nch,'single');
    psdm_evk{con} = nan(nf,nch);
    if nresp < 2
        run_notes(end+1) = string(condis{con}) + ": fewer than 2 response trials, no spectra";
        psdm_tot{con} = nan(nf,nch);
        psdm_ind{con} = nan(nf,nch);
        psdm_pre{con} = nan(nf,nch);
        continue
    end

    for ch = 1:nch
        X = double(reshape(spkcount{con}(rtrs,psdidx,ch), nresp, nwin))'; % [samples x trials] counts, response window
        Xpre = double(reshape(spkcount{con}(rtrs,psdpreidx,ch), nresp, nwin))'; % [samples x trials] counts, pre-stim

        m = mean(X,1); % [1 x trials] mean count per bin of each trial
        mpre = mean(Xpre,1);
        poiss = 2*m/fs; %spectrum level of random (Poisson) firing at each trial's rate (one-sided)
        poisspre = 2*mpre/fs;
        poiss(poiss == 0) = NaN; %a trial without spikes has no spectrum
        poisspre(poisspre == 0) = NaN;

        %total: each trial minus its own mean
        ptot = pmtm(X - m, psd_nw, nfft, fs); % [nfft/2+1 x trials], one column per trial
        psd_tot{con}(:,:,ch) = (ptot(fkeep,:) ./ poiss)';

        %induced: each trial minus the mean of the other trials (leave-one-out, so a trial doesn't take out part of itself)
        loo = (sum(X,2) - X) / (nresp-1); % [samples x trials] mean of the other trials
        Xind = X - loo;
        pind = pmtm(Xind - mean(Xind,1), psd_nw, nfft, fs);
        psd_ind{con}(:,:,ch) = (pind(fkeep,:) ./ poiss)';

        %pre-stim
        ppre = pmtm(Xpre - mpre, psd_nw, nfft, fs);
        psd_pre{con}(:,:,ch) = (ppre(fkeep,:) ./ poisspre)';

        %evoked: the condition mean trace, on the same scale as total
        mtr = mean(X,2); % [samples x 1]
        pevk = pmtm(mtr - mean(mtr), psd_nw, nfft, fs);
        psdm_evk{con}(:,ch) = pevk(fkeep) / (2*mean(m)/fs);
    end

    %line noise out (if asked)
    psd_tot{con}(:,linemask,:) = NaN;
    psd_ind{con}(:,linemask,:) = NaN;
    psd_pre{con}(:,linemask,:) = NaN;
    psdm_evk{con}(linemask,:) = NaN;

    %condition means of the trial spectra
    psdm_tot{con} = reshape(mean(psd_tot{con},1,'omitnan'), nf, nch);
    psdm_ind{con} = reshape(mean(psd_ind{con},1,'omitnan'), nf, nch);
    psdm_pre{con} = reshape(mean(psd_pre{con},1,'omitnan'), nf, nch);
end
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data Visualization: Spike-rate spectra in the probe layout

%one figure per spectrum (total, pre-stim, evoked, induced), every channel in the probe layout, the
%conditions overlaid (group 1 blues, group 2 reds). Dashed line at 1 = random firing, bumps above it = rhythm
if plot_psd == 1
    disp('MUA: Making probe-layout spectrum pictures')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    [nrow,nshank] = size(pmap);

    %colour of each condition: group 1 blues, group 2 reds, subgroups go from dark to light
    grpdark = [0 0.25 0.6; 0.6 0.1 0.1];
    grplight = [0.45 0.7 1; 1 0.55 0.45];
    condcol = zeros(ncond,3);
    for g = 1:numel(condis_grp)
        gcons = find(condgrp == g); %this group's conditions
        for i = 1:numel(gcons)
            if numel(gcons) == 1
                s = 0;
            else
                s = (i-1)/(numel(gcons)-1); %0 = darkest, 1 = lightest
            end
            condcol(gcons(i),:) = grpdark(g,:)*(1-s) + grplight(g,:)*s;
        end
    end

    %legend names (condition + how many response trials went in)
    legnames = cell(1,ncond);
    for con = 1:ncond
        legnames{con} = [condis{con} ' (' num2str(sum(resp{con})) ' trials)'];
    end

    %the four spectra, what they mean, and their file name tags
    specs = {psdm_tot, psdm_pre, psdm_evk, psdm_ind};
    snames = {['total: all rhythm in the response (' num2str(psd_ms(1)) '-' num2str(psd_ms(2)) ' ms)'], ...
        ['pre-stim: rhythm before the stimulus (' num2str(-diff(psd_ms)) ' to 0 ms)'], ...
        'evoked: rhythm locked to the stimulus (spectrum of the mean trace; random-firing floor ~1/trials)', ...
        'induced: rhythm NOT locked to the stimulus (each trial minus the mean of the other trials)'};
    stags = {'tot','pre','evk','ind'};

    for s = 1:numel(specs)
        %one y range for the whole probe
        ytop = 0;
        for con = 1:ncond
            ytop = max(ytop, max(specs{s}{con},[],'all','omitnan'));
        end
        if isnan(ytop) || ytop == 0
            ytop = 1;
        end
        ytop = ytop*1.05;

        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout(nrow,nshank,'TileSpacing','none','Padding','tight');
        for row = 1:nrow
            for shank = 1:nshank
                chid = pmap(row,shank); %raw channel ID at this spot
                ax = nexttile(tl,(row-1)*nshank + shank);
                if chid == 0 %no site here
                    axis(ax,'off')
                    continue
                end
                ch = find(chanlist == chid);
                if isempty(ch) %removed channel
                    axis(ax,'off')
                    text(ax,0.5,0.5,['CH' num2str(chid) ' removed'],'Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5])
                    continue
                end
                axlast = ax; %last panel with data, for the legend
                hold(ax,'on')
                for con = 1:ncond
                    plot(ax, psd_f, specs{s}{con}(:,ch), 'Color',condcol(con,:), 'LineWidth',1)
                end
                if s ~= 3 %random firing = 1 (the evoked floor is ~1/trials instead)
                    yline(ax,1,'k--')
                end
                hold(ax,'off')
                xlim(ax,psd_hz)
                ylim(ax,[0 ytop])
                if respchans(ch) == true %responsive channel: bold blue label
                    text(ax,0.98,0.9,['CH' num2str(chid)],'Units','normalized','HorizontalAlignment','right','Color','b','FontWeight','bold')
                else
                    text(ax,0.98,0.9,['CH' num2str(chid)],'Units','normalized','HorizontalAlignment','right','Color',[.5 .5 .5])
                end
                if row < nrow
                    set(ax,'XTickLabel',[])
                end
                if shank > 1
                    set(ax,'YTickLabel',[])
                end
            end
        end
        lgd = legend(axlast, legnames, 'Orientation','horizontal', 'Interpreter','none');
        lgd.Layout.Tile = 'south';
        xlabel(tl,'Hz')
        ylabel(tl,'power / random firing (1 = no rhythm)')
        if s ~= 3
            linetxt = ', dashed line = random firing (1)';
        else
            linetxt = '';
        end
        title(tl,{[animal ' ' ProbeInfo.Areas{prb} '-Probe spike-rate spectrum, ' snames{s}], ...
            ['line colours: blues = ' grpnames{1} ' (' strjoin(condis_grp{1},', ') '), reds = ' grpnames{2} ' (' strjoin(condis_grp{2},', ') ')' ...
            linetxt ', bold blue CH = responsive channel']},'Interpreter','none')

        if plotsave == 1
            fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'psd_' stags{s}]);
            savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
        end
    end
end




%Now beginning the trial alignment section

%% Data Analysis: Trial alignment by cross-correlation

%Maybe the trials share a rhythm but start a little earlier or later. Each response trial gets ONE shift
%(lag, the same for every channel), found by cross-correlating its spike rate (summed over the responsive
%channels, in align_ms) with a template, biggest shift max_lag_ms either way. Three templates:
%  loo       - the mean of the OTHER trials, redone woody_iter times as the trials get shifted (Woody's method)
%  strongest - the trial with the highest response score
%  first     - the first response trial of the recording
%Null: the exact same thing on a same-length stretch before the stimulus, where nothing should line up.
%Any alignment makes trials look a bit more alike, even pure noise, so the later sections compare to this.
%Lag convention: a lag of +d ms means the trial ran d ms late, so the aligned trial(t) = trial(t + d)
if do_align == 1
    disp('MUA: Lining the trials up')
    tic

    align_methods = {'loo','strongest','first'};
    nmeth = numel(align_methods);
    maxlag = round(max_lag_ms/bin_ms); %biggest shift, in bins
    %a trial shifted by up to max_lag_ms still has to be inside the loaded data (load_ms), so trim the window
    %by max_lag_ms where it would run off the end (e.g. 350-3000 -> 350-2950 with data up to 3000 ms)
    align_use_ms = [max(align_ms(1), load_ms(1) + max_lag_ms), min(align_ms(2), load_ms(2) - max_lag_ms)]; %ms, window actually used
    if isequal(align_use_ms, align_ms) == false
        run_notes(end+1) = "Alignment window trimmed to fit the loaded data: " + align_use_ms(1) + " to " + align_use_ms(2) + " ms (align_ms was " + align_ms(1) + " to " + align_ms(2) + ")";
    end
    alignpre_ms = [-diff(align_use_ms) 0]; %ms, same-length stretch before the stimulus (null)
    if alignpre_ms(1) - max_lag_ms < load_ms(1)
        error(['The pre-stim copy of the alignment window (' num2str(alignpre_ms(1)) ' to 0 ms, +/- ' num2str(max_lag_ms) ' ms of shift) runs off the loaded data (from ' num2str(load_ms(1)) ' ms). Shorten align_ms'])
    end
    winbins = {find(t_bin >= align_use_ms(1) & t_bin < align_use_ms(2)), find(t_bin >= alignpre_ms(1) & t_bin < alignpre_ms(2))}; %bins of {real window, null window}
    nwb = numel(winbins{1});

    lags = cell(1,nmeth); %{methods}{conditions} [response trials x 1] ms, real window
    lags_null = cell(1,nmeth); %{methods}{conditions} [response trials x 1] ms, pre-stim window
    lagpk = cell(1,nmeth); %{methods}{conditions} [response trials x 1] peak cross-correlation (r) at that lag
    lagsd = nan(nmeth,ncond); %ms, SD of the lags (jitter), real window
    lagsd_null = nan(nmeth,ncond); %ms, the same on the pre-stim window
    for m = 1:nmeth
        lags{m} = cell(1,ncond);
        lags_null{m} = cell(1,ncond);
        lagpk{m} = cell(1,ncond);
    end

    for con = 1:ncond
        rtrs = find(resp{con} == true); %response trials
        nresp = numel(rtrs);
        if nresp < 2
            run_notes(end+1) = string(condis{con}) + ": fewer than 2 response trials, no alignment (lags 0)";
            for m = 1:nmeth
                lags{m}{con} = zeros(nresp,1);
                lags_null{m}{con} = zeros(nresp,1);
                lagpk{m}{con} = nan(nresp,1);
            end
            continue
        end
        sig = double(sum(rate{con}(rtrs,:,respchans),3)); % [response trials x all bins] spikes/s summed over the responsive channels
        [~,strongest] = max(respscore{con}(rtrs)); %response trial with the highest score

        for w = 1:2 %1 = real window, 2 = null (pre-stim) window
            wb = winbins{w}; %bins of this window
            for m = 1:nmeth
                trlag = zeros(nresp,1); %bins, current lag of each trial
                trpk = nan(nresp,1);
                if m == 1
                    nround = woody_iter; %leave-one-out: redo as the trials get shifted
                else
                    nround = 1; %fixed template: once is enough
                end
                for it = 1:nround
                    %template for this round
                    if m == 1
                        shifted = zeros(nresp,nwb); %every trial's window at its current shift
                        for tr = 1:nresp
                            shifted(tr,:) = sig(tr, wb + trlag(tr));
                        end
                        shiftsum = sum(shifted,1);
                    elseif m == 2
                        tref = strongest; %template trial
                        tmpl = sig(tref, wb);
                    else
                        tref = 1; %first response trial
                        tmpl = sig(tref, wb);
                    end

                    newlag = zeros(nresp,1);
                    for tr = 1:nresp
                        if m == 1
                            tmpl = (shiftsum - shifted(tr,:)) / (nresp-1); %mean of the OTHER trials
                        elseif tr == tref %the template trial stays where it is
                            newlag(tr) = 0;
                            trpk(tr) = 1;
                            continue
                        end
                        x = sig(tr,wb) - mean(sig(tr,wb)); %this trial, unshifted
                        y = tmpl - mean(tmpl);
                        [r,lagv] = xcorr(x, y, maxlag, 'coeff'); %r peaks at +d when the trial is d late
                        if all(isnan(r)) %a flat trial or template (no spikes), can't be aligned
                            newlag(tr) = 0;
                            trpk(tr) = NaN;
                            continue
                        end
                        [pk,ip] = max(r);
                        newlag(tr) = lagv(ip);
                        trpk(tr) = pk;
                    end
                    trlag = newlag;
                end

                %store, in ms
                if w == 1
                    lags{m}{con} = trlag * bin_ms;
                    lagpk{m}{con} = trpk;
                    lagsd(m,con) = std(trlag * bin_ms);
                else
                    lags_null{m}{con} = trlag * bin_ms;
                    lagsd_null(m,con) = std(trlag * bin_ms);
                end
            end
        end
    end
    clear sig shifted shiftsum
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: Compare the alignment methods

%picture 1: how far each trial got shifted, per method and condition, next to the same thing on the
%pre-stim null. Real structure = lags bunched up (small SD), noise = lags spread out toward +/- max_lag_ms.
%picture 2 (per condition): the responsive-channel summed rate of every response trial, unaligned and
%aligned by each method, with the mean on top. A good alignment makes the mean sharper / taller.
if do_align == 1 && plot_align == 1
    disp('MUA: Making alignment comparison pictures')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    grpdark = [0 0.25 0.6; 0.6 0.1 0.1]; %group 1 blue, group 2 red
    edges = -max_lag_ms-bin_ms/2 : 5 : max_lag_ms+bin_ms/2; %ms, histogram bins (5 ms wide)

    %picture 1: lag histograms, rows = methods, columns = conditions
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(nmeth,ncond,'TileSpacing','compact','Padding','compact');
    for m = 1:nmeth
        for con = 1:ncond
            ax = nexttile(tl,(m-1)*ncond + con);
            histogram(ax, lags_null{m}{con}, edges, 'FaceColor',[.7 .7 .7], 'EdgeColor','none') %null (pre-stim)
            hold(ax,'on')
            histogram(ax, lags{m}{con}, edges, 'FaceColor',grpdark(condgrp(con),:), 'FaceAlpha',.6, 'EdgeColor','none') %real
            xline(ax,0,'k')
            hold(ax,'off')
            xlim(ax,[-max_lag_ms max_lag_ms])
            title(ax,[align_methods{m} ', ' condis{con} ': jitter SD ' num2str(lagsd(m,con),'%.1f') ' ms (null ' num2str(lagsd_null(m,con),'%.1f') ')'],'Interpreter','none')
            if m == nmeth
                xlabel(ax,'lag (ms, + = trial ran late)')
            end
            if con == 1
                ylabel(ax,'trials')
            end
        end
    end
    lgd = legend(ax,{'null (pre-stim)','real'},'Orientation','horizontal');
    lgd.Layout.Tile = 'south';
    title(tl,[animal ' ' ProbeInfo.Areas{prb} '-Probe: trial lags per alignment method (window ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms, max +/- ' num2str(max_lag_ms) ' ms)'],'Interpreter','none')

    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'align_lags']);
        savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
    end

    %picture 2: the trials unaligned and aligned by each method, per condition
    pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
    t_plot = t_bin(pidx); %ms
    pidx_ok = pidx(pidx > maxlag & pidx <= numel(t_bin) - maxlag); %bins that stay inside the data when shifted
    t_ok = t_bin(pidx_ok);
    for con = 1:ncond
        rtrs = find(resp{con} == true); %response trials
        nresp = numel(rtrs);
        if nresp < 2
            continue
        end
        sig = double(sum(rate{con}(rtrs,:,respchans),3)); % [response trials x all bins] spikes/s summed over the responsive channels

        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout(1,nmeth+1,'TileSpacing','compact','Padding','compact');
        panelnames = [{'unaligned'} align_methods];
        ytop = max(sig(:,pidx),[],'all') * 1.05;
        for p = 1:nmeth+1
            if p == 1
                trlag = zeros(nresp,1); %bins
            else
                trlag = round(lags{p-1}{con} / bin_ms); %bins
            end
            sh = zeros(nresp,numel(pidx_ok)); % [response trials x bins shown] each trial read at t + its lag
            for tr = 1:nresp
                sh(tr,:) = sig(tr, pidx_ok + trlag(tr));
            end
            shmean = mean(sh,1);
            aidx = t_ok >= align_use_ms(1) & t_ok < align_use_ms(2); %alignment window within the shown bins

            ax = nexttile(tl,p);
            plot(ax, t_ok, sh', 'Color',[.7 .7 .7], 'LineWidth',0.5) %every trial
            hold(ax,'on')
            plot(ax, t_ok, shmean, 'k', 'LineWidth',1.5) %their mean
            xline(ax,0,'r')
            xline(ax,align_use_ms(1),'b--')
            hold(ax,'off')
            xlim(ax,plot_ms)
            ylim(ax,[0 ytop])
            xlabel(ax,'ms from stim onset')
            if p == 1
                ylabel(ax,'spikes/s, summed over responsive channels')
            end
            %how wiggly the mean is in the alignment window (SD over time): goes up when trials line up better
            title(ax,{panelnames{p}, ['mean''s SD in the window: ' num2str(std(shmean(aidx)),'%.1f') ' spikes/s']},'Interpreter','none')
        end
        title(tl,[animal ' ' ProbeInfo.Areas{prb} '-Probe ' condis{con} ': response trials (grey) and mean (black), before and after alignment (blue dashed = start of the alignment window)'],'Interpreter','none')

        if plotsave == 1
            fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-' 'align_trials']);
            savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
        end
    end
    clear sig sh
end




%Now beginning the aligned-trial analyses

%% Data Analysis: Spectra and trial-to-trial correlation on aligned trials

%Redo the spectra and the trial consistency with the trials lined up (strongest-trial template from the
%alignment section), next to the same thing unaligned, in the alignment window (align_use_ms).
%Spectra: only evoked (spectrum of the condition mean) and induced (trial minus the mean of the OTHER
%trials) are redone. A trial's own (total) spectrum doesn't change when the trial is shifted in time.
%Trial-to-trial correlation (as NNx_LFPCorrelation_LvsLW.m): every pair of response trials, their smoothed
%rate time courses correlated (r), Fisher z, averaged within each condition and within each group
%(subgroups pooled), back to r. Group 2 vs group 1: the group labels are shuffled n_perm times over the
%same trial pairs, p = how often a shuffle gives a difference at least as big as the real one.
%Null: the same consistency numbers on the same-length pre-stim stretch, aligned with the null lags, so
%you can see how much "consistency" turns up with no stimulus (lining trials up always adds some).
if do_align == 1
    disp('MUA: Spectra and trial consistency, unaligned vs aligned')
    tic

    al_names = {'unaligned','strongest'}; %the versions compared
    al_meth = [0 find(strcmp(align_methods,'strongest'))]; %which alignment method's lags each version uses (0 = no shift)
    nver = numel(al_names);
    fs = 1000/bin_ms; %Hz, sampling rate of the counts
    wbr = winbins{1}; %bins of the alignment window (from the alignment section)
    wbn = winbins{2}; %bins of its pre-stim copy (null)
    nwin_al = numel(wbr);

    %spectrum frequencies for this window length
    nfft_al = 2^nextpow2(nwin_al);
    f_al = (0:nfft_al/2)' * fs/nfft_al; %Hz
    fkeep_al = f_al >= psd_hz(1) & f_al <= psd_hz(2);
    psdal_f = f_al(fkeep_al); % [freqs x 1] Hz
    nf_al = numel(psdal_f);
    linemask_al = false(nf_al,1);
    if line_out == 1
        for b = 1:size(line_hz,1)
            linemask_al(psdal_f >= line_hz(b,1) & psdal_f <= line_hz(b,2)) = true;
        end
    end

    psdal_evk = cell(1,nver); %{versions}{conditions} [freqs x channels] evoked spectrum, normalized
    psdal_ind = cell(1,nver); %{versions}{conditions} [freqs x channels] induced spectrum (mean over trials), normalized
    tc_cond = cell(1,nver); %{versions} [channels x conditions] mean pairwise trial r within each condition
    tc_trial = cell(1,nver); %{versions}{conditions} [response trials x channels] each trial's mean r with the other trials of its condition
    tc_grp = cell(1,nver); %{versions} [channels x groups] mean pairwise trial r within each group
    tc_p = cell(1,nver); %{versions} [channels x 1] permutation p, group 2 vs group 1
    tc_cond_null = cell(1,nver); %{versions} [channels x conditions] the same on the pre-stim window
    tc_grp_null = cell(1,nver); %{versions} [channels x groups]
    ngrp = numel(condis_grp);

    for v = 1:nver
        disp(['MUA: ' al_names{v}])
        psdal_evk{v} = cell(1,ncond);
        psdal_ind{v} = cell(1,ncond);
        tc_trial{v} = cell(1,ncond);

        %line every condition's response trials up (or not) and stack them
        allr = []; % [all response trials x window bins x channels] smoothed rate, real window
        allrn = []; % [all response trials x window bins x channels] smoothed rate, pre-stim window
        condlab = []; %condition of each stacked trial
        for con = 1:ncond
            rtrs = find(resp{con} == true); %response trials
            nresp = numel(rtrs);
            psdal_evk{v}{con} = nan(nf_al,nch);
            psdal_ind{v}{con} = nan(nf_al,nch);
            tc_trial{v}{con} = nan(nresp,nch);
            if nresp < 2
                continue
            end
            if al_meth(v) == 0
                lb = zeros(nresp,1); %no shift
                lbn = zeros(nresp,1);
            else
                lb = round(lags{al_meth(v)}{con} / bin_ms); %bins
                lbn = round(lags_null{al_meth(v)}{con} / bin_ms);
            end
            r_al = zeros(nresp,nwin_al,nch,'single'); % [response trials x window bins x channels] rate, each trial read at t + its lag
            c_al = zeros(nresp,nwin_al,nch,'single'); %the same for the raw counts
            r_aln = zeros(nresp,nwin_al,nch,'single'); %rate, pre-stim window
            for tr = 1:nresp
                r_al(tr,:,:) = rate{con}(rtrs(tr), wbr + lb(tr), :);
                c_al(tr,:,:) = spkcount{con}(rtrs(tr), wbr + lb(tr), :);
                r_aln(tr,:,:) = rate{con}(rtrs(tr), wbn + lbn(tr), :);
            end
            allr = cat(1, allr, r_al);
            allrn = cat(1, allrn, r_aln);
            condlab = [condlab; con*ones(nresp,1)];

            %spectra: evoked and induced (raw counts, same normalization as the spectrum section)
            for ch = 1:nch
                X = double(reshape(c_al(:,:,ch), nresp, nwin_al))'; % [samples x trials] counts
                m = mean(X,1); % [1 x trials] mean count per bin
                poiss = 2*m/fs; %random-firing level of each trial
                poiss(poiss == 0) = NaN;
                loo = (sum(X,2) - X) / (nresp-1); %mean of the other trials
                Xind = X - loo;
                pind = pmtm(Xind - mean(Xind,1), psd_nw, nfft_al, fs);
                psdal_ind{v}{con}(:,ch) = mean(pind(fkeep_al,:) ./ poiss, 2, 'omitnan');
                mtr = mean(X,2); %condition mean trace
                pevk = pmtm(mtr - mean(mtr), psd_nw, nfft_al, fs);
                psdal_evk{v}{con}(:,ch) = pevk(fkeep_al) / (2*mean(m)/fs);
            end
            psdal_evk{v}{con}(linemask_al,:) = NaN;
            psdal_ind{v}{con}(linemask_al,:) = NaN;
        end
        grplab = condgrp(condlab)'; %group of each stacked trial
        ntot = numel(condlab);

        %trial-to-trial correlation, per channel
        tc_cond{v} = nan(nch,ncond);
        tc_grp{v} = nan(nch,ngrp);
        tc_p{v} = nan(nch,1);
        tc_cond_null{v} = nan(nch,ncond);
        tc_grp_null{v} = nan(nch,ngrp);
        rng('default') %same shuffles every run
        for ch = 1:nch
            for w = 1:2 %1 = real window, 2 = pre-stim (null) window
                if w == 1
                    X = double(reshape(allr(:,:,ch), ntot, nwin_al)); % [trials x bins]
                else
                    X = double(reshape(allrn(:,:,ch), ntot, nwin_al));
                end
                C = corr(X'); % [trials x trials] r of every pair of trial time courses
                C(C > 0.9999) = 0.9999; %so atanh stays finite
                Z = atanh(C); %Fisher z
                Z(1:ntot+1:end) = NaN; %leave out each trial with itself

                %per condition
                for con = 1:ncond
                    idx = condlab == con;
                    if sum(idx) < 2
                        continue
                    end
                    if w == 1
                        tc_cond{v}(ch,con) = tanh(mean(Z(idx,idx),'all','omitnan'));
                        tc_trial{v}{con}(:,ch) = tanh(mean(Z(idx,idx),2,'omitnan')); %each trial with the others of its condition
                    else
                        tc_cond_null{v}(ch,con) = tanh(mean(Z(idx,idx),'all','omitnan'));
                    end
                end

                %per group
                zg = nan(1,ngrp);
                for g = 1:ngrp
                    gi = grplab == g;
                    if sum(gi) >= 2
                        zg(g) = mean(Z(gi,gi),'all','omitnan');
                    end
                end
                if w == 1
                    tc_grp{v}(ch,:) = tanh(zg);
                else
                    tc_grp_null{v}(ch,:) = tanh(zg);
                end

                %group 2 vs group 1, label shuffle (real window only)
                if w == 1 && ngrp == 2 && all(isnan(zg) == false)
                    obs = zg(2) - zg(1);
                    nulldiff = zeros(n_perm,1);
                    for p = 1:n_perm
                        shlab = grplab(randperm(ntot)); %shuffled group labels
                        g1 = shlab == 1;
                        g2 = shlab == 2;
                        nulldiff(p) = mean(Z(g2,g2),'all','omitnan') - mean(Z(g1,g1),'all','omitnan');
                    end
                    tc_p{v}(ch) = (sum(abs(nulldiff) >= abs(obs)) + 1) / (n_perm + 1);
                end
            end
        end
    end
    clear allr allrn r_al c_al r_aln X C Z
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: Aligned-trial results

%picture 1: trial consistency on the probe (one square per channel): group 1 r, group 2 r, group 2 - 1
%           (black square = p < 0.05), and the same difference on the pre-stim null; unaligned and aligned
%picture 2: how much lining the trials up raised the consistency, per group, next to the same gain on the null
%picture 3: evoked and induced spectra, unaligned (dashed) vs aligned (solid), probe layout
%picture 4: each trial's consistency with its condition (mean over the responsive channels), unaligned vs aligned
if do_align == 1 && (plot_corr == 1 || plot_psd == 1)
    disp('MUA: Making aligned-trial pictures')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    [nrow,nshank] = size(pmap);
    vun = find(strcmp(al_names,'unaligned'));
    vstr = find(strcmp(al_names,'strongest'));
    [rrow,rshank] = find(pmap > 0 & ismember(pmap, chanlist) == false); %removed channels
    nc = 128; %blue - white - red colour map for differences (0 = white)
    divmap = [linspace(0,1,nc)' linspace(0.3,1,nc)' ones(nc,1); ones(nc,1) linspace(1,0.2,nc)' linspace(1,0.2,nc)'];

    %colour of each condition: group 1 blues, group 2 reds, subgroups dark to light (as the spectrum pictures)
    grpdark = [0 0.25 0.6; 0.6 0.1 0.1];
    grplight = [0.45 0.7 1; 1 0.55 0.45];
    condcol = zeros(ncond,3);
    for g = 1:numel(condis_grp)
        gcons = find(condgrp == g);
        for i = 1:numel(gcons)
            if numel(gcons) == 1
                s = 0;
            else
                s = (i-1)/(numel(gcons)-1);
            end
            condcol(gcons(i),:) = grpdark(g,:)*(1-s) + grplight(g,:)*s;
        end
    end

    if plot_corr == 1
        %pictures 1 and 2: probe maps. Each figure is a list of panels: values (channels x 1), title,
        %diverging colours or not, channels to mark with a black square
        for fg = 1:2
            pv = {};
            pt = {};
            pd = [];
            pm = {};
            if fg == 1
                for v = 1:nver
                    pv(end+1:end+4) = {tc_grp{v}(:,1), tc_grp{v}(:,2), tc_grp{v}(:,2) - tc_grp{v}(:,1), tc_grp_null{v}(:,2) - tc_grp_null{v}(:,1)};
                    pt(end+1:end+4) = {[al_names{v} ': ' grpnames{1} ' trial r'], [al_names{v} ': ' grpnames{2} ' trial r'], ...
                        [al_names{v} ': ' grpnames{2} ' - ' grpnames{1} ' (square = p < 0.05)'], [al_names{v} ': ' grpnames{2} ' - ' grpnames{1} ', pre-stim null']};
                    pd(end+1:end+4) = [0 0 1 1];
                    pm(end+1:end+4) = {[], [], tc_p{v} < 0.05, []};
                end
                tlsize = [nver 4];
                figtitle = [animal ' ' ProbeInfo.Areas{prb} '-Probe: trial-to-trial consistency (mean pairwise r of the smoothed rate, ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms)'];
                ftag = 'tc_maps';
            else
                for g = 1:2
                    pv{end+1} = tc_grp{vstr}(:,g) - tc_grp{vun}(:,g);
                    pt{end+1} = [grpnames{g} ': gain from aligning (aligned - unaligned r)'];
                end
                for g = 1:2
                    pv{end+1} = tc_grp_null{vstr}(:,g) - tc_grp_null{vun}(:,g);
                    pt{end+1} = [grpnames{g} ': same gain on the pre-stim null'];
                end
                pd = [1 1 1 1];
                pm = {[], [], [], []};
                tlsize = [1 4];
                figtitle = [animal ' ' ProbeInfo.Areas{prb} '-Probe: how much lining the trials up raised their consistency (real vs no-stim null)'];
                ftag = 'tc_gain';
            end

            %colour ranges: one for the r panels, one symmetric one for the difference panels
            rvals = cell2mat(pv(pd == 0)');
            dvals = cell2mat(pv(pd == 1)');
            rlim = [min(rvals,[],'omitnan') max(rvals,[],'omitnan')];
            if isempty(rlim) || any(isnan(rlim)) || rlim(1) == rlim(2)
                rlim = [0 1];
            end
            dmax = max(abs(dvals),[],'omitnan');
            if isempty(dmax) || isnan(dmax) || dmax == 0
                dmax = 0.01;
            end

            figure()
            set(gcf, 'Position', get(0, 'Screensize'));
            tl = tiledlayout(tlsize(1),tlsize(2),'TileSpacing','compact','Padding','compact');
            for k = 1:numel(pv)
                ax = nexttile(tl,k);
                zmap = nan(size(pmap)); % [rows x shanks]
                for ch = 1:nch
                    zmap(chrow(ch),chshank(ch)) = pv{k}(ch);
                end
                h = imagesc(ax, zmap);
                set(h,'AlphaData',isnan(zmap) == false) %no site / removed = blank
                if pd(k) == 1
                    colormap(ax,divmap)
                    clim(ax,[-dmax dmax])
                else
                    colormap(ax,parula)
                    clim(ax,rlim)
                end
                colorbar(ax)
                hold(ax,'on')
                for ch = 1:nch
                    if respchans(ch) == true
                        text(ax, chshank(ch), chrow(ch), num2str(chanlist(ch)), 'HorizontalAlignment','center','FontSize',7,'FontWeight','bold')
                    else
                        text(ax, chshank(ch), chrow(ch), num2str(chanlist(ch)), 'HorizontalAlignment','center','FontSize',7,'Color',[.4 .4 .4])
                    end
                end
                if isempty(pm{k}) == false && any(pm{k})
                    plot(ax, chshank(pm{k}), chrow(pm{k}), 'ks', 'MarkerSize',16, 'LineWidth',1.5) %p < 0.05
                end
                for i = 1:numel(rrow)
                    text(ax, rshank(i), rrow(i), 'x', 'HorizontalAlignment','center','Color',[.5 .5 .5])
                end
                hold(ax,'off')
                set(ax,'XTick',1:nshank,'YTick',1:nrow)
                xlabel(ax,'shank')
                ylabel(ax,'row (1 = top)')
                title(ax, pt{k}, 'Interpreter','none')
            end
            title(tl,{figtitle, 'bold = responsive channel, x = removed'},'Interpreter','none')

            if plotsave == 1
                fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' ftag]);
                savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
            end
        end

        %picture 4: each trial's consistency with its condition, averaged over the responsive channels
        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        ax = axes();
        hold(ax,'on')
        xt = [];
        xtl = {};
        for con = 1:ncond
            trval = cell(1,nver); %each trial's value per version
            for v = 1:nver
                trval{v} = mean(tc_trial{v}{con}(:,respchans),2,'omitnan'); % [response trials x 1]
            end
            xpos = (con-1)*3 + [1 2]; %unaligned, aligned
            ntr = numel(trval{1});
            if ntr == 0
                continue
            end
            %lines joining each trial's unaligned and aligned value
            plot(ax, xpos', [trval{vun} trval{vstr}]', 'Color',[.75 .75 .75], 'LineWidth',0.5)
            for v = 1:nver
                if v == vun
                    col = condcol(con,:)*0.4 + [1 1 1]*0.6; %lighter for unaligned
                else
                    col = condcol(con,:);
                end
                boxchart(ax, xpos(v)*ones(ntr,1), trval{v}, 'BoxFaceColor',col, 'MarkerStyle','none', 'BoxWidth',0.5)
                scatter(ax, xpos(v) + (rand(ntr,1)-0.5)*0.25, trval{v}, 12, col, 'filled') %each trial
            end
            xt = [xt xpos];
            xtl = [xtl {[condis{con} ' unaligned'], [condis{con} ' aligned']}];
        end
        hold(ax,'off')
        set(ax,'XTick',xt,'XTickLabel',xtl,'TickLabelInterpreter','none')
        ylabel(ax,'trial''s mean r with the other trials of its condition (mean over responsive channels)')
        title(ax,{[animal ' ' ProbeInfo.Areas{prb} '-Probe: per-trial consistency, unaligned vs aligned (strongest-trial template)'], ...
            'grey lines join the same trial before and after aligning'},'Interpreter','none')

        if plotsave == 1
            fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'tc_trials']);
            savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
        end
    end

    if plot_psd == 1
        %picture 3: evoked and induced spectra, unaligned (dashed) vs aligned (solid), probe layout
        specs = {psdal_evk, psdal_ind};
        snames = {'evoked (spectrum of the condition mean; random-firing floor ~1/trials)', 'induced (each trial minus the mean of the other trials)'};
        stags = {'psdal_evk','psdal_ind'};
        lstyle = {'--','-'}; %unaligned dashed, aligned solid (order of al_names)
        for s = 1:2
            ytop = 0;
            for v = 1:nver
                for con = 1:ncond
                    ytop = max(ytop, max(specs{s}{v}{con},[],'all','omitnan'));
                end
            end
            if isnan(ytop) || ytop == 0
                ytop = 1;
            end
            ytop = ytop*1.05;

            figure()
            set(gcf, 'Position', get(0, 'Screensize'));
            tl = tiledlayout(nrow,nshank,'TileSpacing','none','Padding','tight');
            for row = 1:nrow
                for shank = 1:nshank
                    chid = pmap(row,shank); %raw channel ID at this spot
                    ax = nexttile(tl,(row-1)*nshank + shank);
                    if chid == 0 %no site here
                        axis(ax,'off')
                        continue
                    end
                    ch = find(chanlist == chid);
                    if isempty(ch) %removed channel
                        axis(ax,'off')
                        text(ax,0.5,0.5,['CH' num2str(chid) ' removed'],'Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5])
                        continue
                    end
                    axlast = ax; %last panel with data, for the legend
                    hold(ax,'on')
                    legnames = {};
                    for con = 1:ncond
                        for v = 1:nver
                            plot(ax, psdal_f, specs{s}{v}{con}(:,ch), lstyle{v}, 'Color',condcol(con,:), 'LineWidth',1)
                            legnames{end+1} = [condis{con} ' ' al_names{v}];
                        end
                    end
                    if s == 2 %random firing = 1 (not for evoked)
                        yline(ax,1,'k:')
                    end
                    hold(ax,'off')
                    xlim(ax,psd_hz)
                    ylim(ax,[0 ytop])
                    if respchans(ch) == true
                        text(ax,0.98,0.9,['CH' num2str(chid)],'Units','normalized','HorizontalAlignment','right','Color','b','FontWeight','bold')
                    else
                        text(ax,0.98,0.9,['CH' num2str(chid)],'Units','normalized','HorizontalAlignment','right','Color',[.5 .5 .5])
                    end
                    if row < nrow
                        set(ax,'XTickLabel',[])
                    end
                    if shank > 1
                        set(ax,'YTickLabel',[])
                    end
                end
            end
            lgd = legend(axlast, legnames, 'Orientation','horizontal', 'Interpreter','none');
            lgd.Layout.Tile = 'south';
            xlabel(tl,'Hz')
            ylabel(tl,'power / random firing')
            title(tl,{[animal ' ' ProbeInfo.Areas{prb} '-Probe ' snames{s} ', ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms'], ...
                ['dashed = unaligned, solid = aligned (strongest-trial template); blues = ' grpnames{1} ', reds = ' grpnames{2} '; bold blue CH = responsive']},'Interpreter','none')

            if plotsave == 1
                fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' stags{s}]);
                savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
            end
        end
    end
end




%Now beginning the S1 section
%(the spike-LFP and V1-S1 LFP phase locking that used to be here were dropped, see commit 1747056)

%% Data Extraction: Load the S1 LFP of the same recordings

%the LFP files sit next to the spiking files (...\LFP\<animal>\<animal>-<stim>-LFP.mat), so no picking.
%Only the S1 probe's kept channels are kept, same response trials in the same order as the MUA.
%Also picks the single S1 channel: kept channel number s1_ref_row (from the top) of the S1 shank with the
%highest s1_shank_label score
if do_S1 == 1
    disp('S1: Loading the S1 LFP of the same recordings')
    tic

    % s1pre_ms = [-diff(psd_ms) 0]; %commented out: the S1 comparison is time-frequency now (sliding windows),
    % s1pre_ms = [-diff(s1_ms) 0]; %its windows are checked against lfp_load_ms in that section
    lfpsamp = stim_onset_lfp + (lfp_load_ms(1):lfp_load_ms(2)-1); %LFP samples to keep (1 kHz)
    t_lfp = lfp_load_ms(1):lfp_load_ms(2)-1; %ms of each kept sample

    %S1 probe channels
    prb_s1 = find(strcmp(string(ProbeInfo.Areas), s1_region));
    if numel(prb_s1) ~= 1
        error(['s1_region "' char(s1_region) '" does not pick exactly one probe in ProbeInfo.Areas'])
    end
    pmap_s1 = ProbeInfo.ProbeMaps{prb_s1}; % [rows x shanks] raw channel IDs
    chremove_s1 = [];
    if numel(ProbeInfo.Ch_Remove) >= prb_s1
        chremove_s1 = ProbeInfo.Ch_Remove{prb_s1};
    end
    chanlist_s1 = setdiff(nonzeros(pmap_s1), chremove_s1(:)); %sorted column of kept S1 raw channel IDs

    %the single S1 channel: kept channel number s1_ref_row of the S1 shank with the highest s1_shank_label score
    if isfield(ProbeInfo,'ShankLabels') == false || numel(ProbeInfo.ShankLabels) < prb_s1 || isstruct(ProbeInfo.ShankLabels{prb_s1}) == false ...
            || isfield(ProbeInfo.ShankLabels{prb_s1}, char(s1_shank_label)) == false
        haslabels = '';
        if isfield(ProbeInfo,'ShankLabels') && numel(ProbeInfo.ShankLabels) >= prb_s1 && isstruct(ProbeInfo.ShankLabels{prb_s1})
            haslabels = strjoin(fieldnames(ProbeInfo.ShankLabels{prb_s1}),', ');
        end
        error(['No ShankLabels.' char(s1_shank_label) ' on the S1 probe. Labels there: ' haslabels])
    end
    shscore = ProbeInfo.ShankLabels{prb_s1}.(char(s1_shank_label)); % [1 x shanks]
    [~,pkshank] = max(shscore); %shank with the highest score
    keptsh = pmap_s1(:, pkshank); %that shank, top to bottom
    keptsh = keptsh(keptsh > 0 & ismember(keptsh, chanlist_s1)); %its kept channels only, still top to bottom
    if numel(keptsh) < s1_ref_row
        error(['S1 shank ' num2str(pkshank) ' has only ' num2str(numel(keptsh)) ' kept channels, s1_ref_row is ' num2str(s1_ref_row)])
    end
    s1_refch = keptsh(s1_ref_row); %raw channel ID
    s1_refidx = find(chanlist_s1 == s1_refch); %its place in chanlist_s1

    lfp_s1 = cell(1,ncond); %{conditions} [response trials x samples x S1 channels] uV, single
    for con = 1:ncond
        for file = 1:numel(condi_files{con})
            spkfile = condi_files{con}{file};
            lfpfile = strrep(strrep(spkfile,'\Spiking\','\LFP\'),'-spiking_results.mat','-LFP.mat');
            if isfile(lfpfile) == false
                error(['No LFP file for ' spkfile ': ' lfpfile])
            end
            [~,fname] = fileparts(lfpfile);
            disp(['S1: Loading ' condis{con} ', ' fname])

            d = load(lfpfile,'stim_lfp_stimchunks','tr_keep','tr_remove'); %(not -v7.3, so the whole thing loads)
            sdat = matfile(spkfile);
            sk = sdat.tr_keep; %read whole first (matfile can't index with (:)), they're small
            sr = sdat.tr_remove;
            if isequal(d.tr_keep(:), sk(:)) == false || isequal(d.tr_remove(:), sr(:)) == false
                error([fname ': tr_keep / tr_remove differ from the spiking file, so the trials would not match'])
            end
            if lfpsamp(1) < 1 || lfpsamp(end) > size(d.stim_lfp_stimchunks,2)
                error([fname ': lfp_load_ms runs off the trials (' num2str(size(d.stim_lfp_stimchunks,2)) ' samples)'])
            end

            trs = trinfo{con}(trinfo{con}(:,1) == file, 2); %this file's good trials, same order as the MUA
            lfp_s1{con} = cat(1, lfp_s1{con}, single(d.stim_lfp_stimchunks(trs, lfpsamp, chanlist_s1)));
            clear d
        end
        lfp_s1{con} = lfp_s1{con}(resp{con},:,:); %keep the response trials (same ones as everything else)
    end
    clear sdat
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Analysis: V1 spike rhythmicity and S1 power over time (sliding windows, group 2)

%The SR idea: if the added input (group 2, e.g. LW) works as "noise" in some band, V1 spiking should get
%more rhythmic in the bands where S1 has extra power, and at the same times. So, for group 2's response
%trials only, in sliding windows (tf_win_ms long, every tf_step_ms, centres tf_ms), per trial:
%  V1: multitaper spectrum of the raw 1 ms spike counts per window, divided by that window's random
%      (Poisson) firing level -> 1 = no rhythm, > 1 = rhythm, averaged over the responsive channels (dB)
%  S1: multitaper LFP power per window for every S1 channel (dB)
%then the shape match of the two time courses at each frequency (below), with a trial-shuffle null
%(trial-averaged maps, r(t) across frequencies and the across-trial r per pixel were tried and dropped)
if do_S1 == 1
    disp(['S1: Sliding-window spectra of V1 spiking and S1 LFP, ' grpnames{2} ' (this takes a minute)'])
    tic

    if bin_ms ~= 1
        error('The V1-S1 time-frequency comparison needs bin_ms = 1 (spike bins = LFP samples, 1 kHz)')
    end
    fs1 = 1000; %Hz, the counts and the LFP
    halfw = tf_win_ms/2; %ms

    %window centres: keep the ones whose windows fit in the loaded MUA (load_ms) and LFP (lfp_load_ms)
    cent = tf_ms(1):tf_step_ms:tf_ms(2);
    fits = cent - halfw >= max(load_ms(1), lfp_load_ms(1)) & cent + halfw <= min(load_ms(2), lfp_load_ms(2));
    if any(fits) == false
        error('No sliding window fits in the loaded data. Change tf_ms / tf_win_ms')
    end
    if all(fits) == false
        run_notes(end+1) = "Time-frequency window centres trimmed to " + min(cent(fits)) + " to " + max(cent(fits)) + " ms (the loaded data end there)";
    end
    tf_t = cent(fits); % [1 x windows] ms, window centres used
    nt = numel(tf_t);

    %multitaper tapers, made once (Slepian sequences, unit energy), and the frequencies kept
    ntap = 2*tf_nw - 1;
    tapers = dpss(tf_win_ms, tf_nw, ntap); % [window samples x tapers]
    nfft_tf = 2^nextpow2(tf_win_ms);
    f_tf_all = (0:nfft_tf/2)' * fs1/nfft_tf; %Hz
    fk_tf = f_tf_all >= psd_hz(1) & f_tf_all <= psd_hz(2);
    tf_f = f_tf_all(fk_tf); % [freqs x 1] Hz
    nftf = numel(tf_f);
    linemask_tf = false(nftf,1);
    if line_out == 1
        for b = 1:size(line_hz,1)
            linemask_tf(tf_f >= line_hz(b,1) & tf_f <= line_hz(b,2)) = true;
        end
    end
    %(one-sided power spectral density of a windowed stretch, averaged over the tapers: same scaling as
    %pmtm, so random firing at mean count m per ms sits at 2*m/fs)

    g2cons = find(condgrp == 2); %group 2's conditions
    rch = find(respchans); %responsive V1 channels
    nrch = numel(rch);
    nch_s1 = numel(chanlist_s1);

    ntr2 = 0; %group 2 response trials

    %per-trial values
    ntr2tot = 0;
    for con = g2cons
        ntr2tot = ntr2tot + size(lfp_s1{con},1);
    end
    v1trsum = zeros(nftf,nt,ntr2tot); % [freqs x windows x trials] V1 rhythmicity summed over the responsive channels
    v1trn = zeros(nftf,nt,ntr2tot); %how many responsive channels went into each
    s1tr = zeros(nftf,nt,ntr2tot,nch_s1,'single'); % [freqs x windows x trials x S1 channels] S1 power, single saves memory
    off = 0; %trials of the earlier conditions

    for con = g2cons
        nresp = size(lfp_s1{con},1); %response trials (lfp_s1 only holds these)
        if nresp < 1
            continue
        end
        rtrs = find(resp{con} == true); %their rows in spkcount
        ntr2 = ntr2 + nresp;
        disp(['S1: ' condis{con}])

        for ti = 1:nt
            b0 = find(t_bin == tf_t(ti) - halfw); %first MUA bin of this window
            l0 = find(t_lfp == tf_t(ti) - halfw); %first LFP sample of this window
            bidx = b0 + (0:tf_win_ms-1);
            lidx = l0 + (0:tf_win_ms-1);

            %V1, responsive channels
            for i = 1:nrch
                X = double(reshape(spkcount{con}(rtrs,bidx,rch(i)), nresp, tf_win_ms))'; % [window samples x trials] spikes per ms
                m = mean(X,1); % [1 x trials] mean count per ms in this window
                P = zeros(nfft_tf/2+1, nresp);
                for j = 1:ntap
                    Y = fft((X - m) .* tapers(:,j), nfft_tf);
                    P = P + abs(Y(1:nfft_tf/2+1,:)).^2;
                end
                P = P / ntap / fs1;
                P(2:end-1,:) = 2*P(2:end-1,:); %one-sided
                lvl = 2*m/fs1; %random-firing level of each trial in this window
                lvl(lvl == 0) = NaN; %no spikes, no rhythmicity
                R = P(fk_tf,:) ./ lvl; % [freqs x trials] 1 = random firing
                good = isnan(R) == false;
                R(good == false) = 0;
                v1trsum(:,ti,off+(1:nresp)) = v1trsum(:,ti,off+(1:nresp)) + reshape(R, nftf, 1, nresp); %per trial
                v1trn(:,ti,off+(1:nresp)) = v1trn(:,ti,off+(1:nresp)) + reshape(good, nftf, 1, nresp);
            end

            %S1, every channel
            for c1 = 1:nch_s1
                X = double(lfp_s1{con}(:,lidx,c1))'; % [window samples x trials] uV
                P = zeros(nfft_tf/2+1, nresp);
                for j = 1:ntap
                    Y = fft((X - mean(X,1)) .* tapers(:,j), nfft_tf);
                    P = P + abs(Y(1:nfft_tf/2+1,:)).^2;
                end
                P = P / ntap / fs1;
                P(2:end-1,:) = 2*P(2:end-1,:);
                s1tr(:,ti,off+(1:nresp),c1) = single(reshape(P(fk_tf,:), nftf, 1, nresp)); %per trial
            end
        end
        off = off + nresp;
    end
    if ntr2 == 0
        error(['No ' grpnames{2} ' response trials, nothing to compare'])
    end

    %V1 rhythmicity per trial: responsive-channel mean, dB above random firing (NaN = no spikes)
    v1tr = 10*log10(v1trsum ./ v1trn); % [freqs x windows x trials]
    v1tr(isinf(v1tr)) = NaN;

    %shape of the time courses: at each frequency, within each trial, correlate the V1 rhythmicity curve and
    %the S1 power curve over a sliding block of windows (shape_block_ms of centres), then average the trials'
    %r (Fisher z). Pearson ignores offset and scale, so a curve sitting below its pre-stim level (negative dB)
    %still counts if it rises and falls with the other one. Neighbouring windows overlap a lot, so each block
    %holds only a few independent points: r is noisy and large values happen by chance -> trial-shuffle null.
    %
    %Trial-shuffle null: V1 curves of each trial paired with ANOTHER trial's S1 curves (no trial keeps its own
    %partner), whole map recomputed, n_shuf times. Whatever survives shuffling (the time course all trials
    %share, the overlapping windows) is chance level. Each pixel gets a z against its shuffles; neighbouring
    %pixels with |z| past the shuf_p threshold form clusters (positive and negative apart), cluster mass =
    %sum of |z|; a real cluster counts when its mass beats the (1 - clust_p) percentile of the biggest chance
    %cluster of each shuffle, over ALL S1 channels (so the whole probe is covered).
    %Speed: each trial's block curves are centred and scaled once (for V1 and for S1), so a pairing is just
    %a multiply and a sum. Windows where V1 had no spikes are left out of the V1 curve but not the S1 curve
    %(slight approximation; rare once the responsive channels are averaged)
    if exist('bwconncomp','file') == 0
        error('The shape-map clusters need bwconncomp (Image Processing Toolbox)')
    end
    hb = round(shape_block_ms/2/tf_step_ms); %windows either side of the block centre
    L = 2*hb + 1; %windows per block
    bc = 1+hb : nt-hb; %block centres (window index)
    nbk = numel(bc);

    %V1: each trial's block curve, centred and scaled to length 1 (once, the same for every S1 channel)
    Vn = zeros(nftf,L,ntr2tot,nbk,'single'); % [freqs x block windows x trials x blocks]
    Vok = false(nftf,1,ntr2tot,nbk); %trial usable in that block (>= 5 windows with spikes, not flat)
    for k = 1:nbk
        a = v1tr(:, bc(k)-hb : bc(k)+hb, :); % [freqs x block windows x trials]
        okv = isnan(a) == false;
        a(okv == false) = 0;
        nv = sum(okv,2);
        da = (a - sum(a,2)./nv) .* okv; %deviations from the trial's block mean (0 where no spikes)
        den = sqrt(sum(da.^2,2));
        Vok(:,:,:,k) = nv >= 5 & den > 0;
        da = da ./ den;
        da(isnan(da)) = 0;
        Vn(:,:,:,k) = single(da);
    end

    %the shuffles: the same pairings for every S1 channel (needed for the max over channels)
    rng('default') %same shuffles every run
    permset = zeros(n_shuf, ntr2tot);
    for s = 1:n_shuf
        pm = randperm(ntr2tot);
        while any(pm == 1:ntr2tot) %no trial paired with itself
            pm = randperm(ntr2tot);
        end
        permset(s,:) = pm;
    end
    zthr = norminv(1 - shuf_p/2); %pixel |z| threshold for the clusters

    r_shape_ch = nan(nftf,nt,nch_s1); % [freqs x windows (block centres) x S1 channels], NaN at the edges
    shape_zpix = nan(nftf,nt,nch_s1); % [freqs x windows x S1 channels] pixel z vs the shuffles
    shape_ppix = nan(nftf,nt,nch_s1); % [freqs x windows x S1 channels] pixel p (two-sided, from the shuffles)
    nullmax = zeros(n_shuf, nch_s1); %biggest chance cluster mass per shuffle and channel
    for c1 = 1:nch_s1
        s1l = 10*log10(double(s1tr(:,:,:,c1))); % [freqs x windows x trials] dB
        s1l(isinf(s1l)) = NaN;

        %S1: each trial's block curve, centred and scaled to length 1
        Sn = zeros(nftf,L,ntr2tot,nbk,'single');
        for k = 1:nbk
            b = s1l(:, bc(k)-hb : bc(k)+hb, :);
            b(isnan(b)) = 0;
            db = b - mean(b,2);
            db = db ./ sqrt(sum(db.^2,2));
            db(isnan(db)) = 0;
            Sn(:,:,:,k) = single(db);
        end

        %real map: each trial's V1 with its own S1, Fisher-z mean over trials
        rt = double(sum(Vn .* Sn, 2)); % [freqs x 1 x trials x blocks] r per trial
        rt(Vok == false) = NaN;
        rt(rt > 0.9999) = 0.9999; %keep atanh finite (NaNs stay NaN)
        rt(rt < -0.9999) = -0.9999;
        zreal = reshape(mean(atanh(rt),3,'omitnan'), nftf, nbk); % [freqs x blocks]

        %shuffled maps
        zsh = zeros(nftf,nbk,n_shuf,'single');
        for s = 1:n_shuf
            rt = double(sum(Vn .* Sn(:,:,permset(s,:),:), 2)); %V1 of trial i with S1 of another trial
            rt(Vok == false) = NaN;
            rt(rt > 0.9999) = 0.9999;
            rt(rt < -0.9999) = -0.9999;
            zsh(:,:,s) = single(reshape(mean(atanh(rt),3,'omitnan'), nftf, nbk));
        end
        zsh(linemask_tf,:,:) = NaN;
        zreal(linemask_tf,:) = NaN;

        %pixel z and p against the shuffles
        mu = mean(zsh,3,'omitnan');
        sd = std(zsh,0,3,'omitnan');
        zp = (zreal - mu) ./ sd; % [freqs x blocks]
        pp = (sum(abs(zsh - mu) >= abs(zreal - mu), 3) + 1) / (n_shuf + 1);
        r_shape_ch(:,bc,c1) = tanh(zreal);
        shape_zpix(:,bc,c1) = zp;
        shape_ppix(:,bc,c1) = pp;

        %biggest chance cluster of each shuffle (its own map against the same shuffle mean / SD)
        for s = 1:n_shuf
            zs = (double(zsh(:,:,s)) - mu) ./ sd;
            mx = 0;
            for sgn = [1 -1] %positive and negative clusters apart
                cc = bwconncomp(sgn*zs > zthr, 4);
                for q = 1:cc.NumObjects
                    mx = max(mx, sum(abs(zs(cc.PixelIdxList{q}))));
                end
            end
            nullmax(s,c1) = mx;
        end
        disp(['S1: shape-map shuffles done for CH' num2str(chanlist_s1(c1)) ' (' num2str(c1) ' of ' num2str(nch_s1) ')'])
    end

    %cluster threshold: biggest chance cluster over all S1 channels, per shuffle
    shape_clustthr = prctile(max(nullmax,[],2), 100*(1 - clust_p));
    shape_sig = false(nftf,nt,nch_s1); %pixels in significant clusters
    for c1 = 1:nch_s1
        zp = shape_zpix(:,bc,c1);
        sigm = false(nftf,nbk);
        for sgn = [1 -1]
            cc = bwconncomp(sgn*zp > zthr, 4);
            for q = 1:cc.NumObjects
                if sum(abs(zp(cc.PixelIdxList{q}))) > shape_clustthr
                    sigm(cc.PixelIdxList{q}) = true;
                end
            end
        end
        shape_sig(:,bc,c1) = sigm;
    end
    clear X P Y R v1trsum v1trn s1tr s1l a b da db rt Vn Sn zsh okv
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: V1 spike rhythmicity vs S1 power over time

%the shape-match maps only: per S1 channel (S1 probe layout), a time-frequency map coloured by how alike the
%shapes of the V1 rhythmicity and S1 power time courses are, with the trial-shuffle significance
%(the V1 / S1 time-frequency maps and the r(t) pictures that were here were dropped as not useful)
if do_S1 == 1 && plot_lfpcoh == 1
    disp('S1: Making the shape-match maps')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    nc = 128; %blue - white - red colour map (0 = white)
    divmap = [linspace(0,1,nc)' linspace(0.3,1,nc)' ones(nc,1); ones(nc,1) linspace(1,0.2,nc)' linspace(1,0.2,nc)'];
    [nrow_s1,nshank_s1] = size(pmap_s1); %S1 probe layout

    %per S1 channel (S1 probe layout), a time-frequency map where the colour = how alike the SHAPE
    %of the V1 rhythmicity and S1 power time courses are at that frequency, around that time (per-trial r
    %over a shape_block_ms block, Fisher-z averaged over group 2's trials)
    rmap = r_shape_ch; %(an across-trial r per pixel was tried first and dropped: it compared values, not the shape over time)
    rlim = max(abs(rmap(:)),[],'omitnan'); %one symmetric colour range for the whole probe
    if isnan(rlim) || rlim == 0
        rlim = 1;
    end
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(nrow_s1,nshank_s1,'TileSpacing','none','Padding','tight');
    for row = 1:nrow_s1
        for shank = 1:nshank_s1
            chid = pmap_s1(row,shank); %raw channel ID at this spot
            ax = nexttile(tl,(row-1)*nshank_s1 + shank);
            if chid == 0 %no site here
                axis(ax,'off')
                continue
            end
            c1 = find(chanlist_s1 == chid);
            if isempty(c1) %removed channel
                axis(ax,'off')
                text(ax,0.5,0.5,['CH' num2str(chid) ' removed'],'Units','normalized','HorizontalAlignment','center','Color',[.5 .5 .5])
                continue
            end
            h = imagesc(ax, tf_t, tf_f, rmap(:,:,c1));
            % set(h,'AlphaData',isnan(rmap(:,:,c1)) == false) %commented out: now faded outside the significant clusters
            set(h,'AlphaData',(isnan(rmap(:,:,c1)) == false) .* (0.35 + 0.65*shape_sig(:,:,c1))) %NaN = blank, not significant = faded
            axis(ax,'xy')
            colormap(ax,divmap)
            clim(ax,[-rlim rlim])
            hold(ax,'on')
            xline(ax,0,'k')
            if any(shape_sig(:,:,c1),'all')
                contour(ax, tf_t, tf_f, double(shape_sig(:,:,c1)), [0.5 0.5], 'k', 'LineWidth',1) %outline the significant clusters
            end
            hold(ax,'off')
            if chid == s1_refch %the single S1 channel
                text(ax,0.02,0.9,['CH' num2str(chid)],'Units','normalized','Color','b','FontWeight','bold')
            else
                text(ax,0.02,0.9,['CH' num2str(chid)],'Units','normalized','Color','k')
            end
            if row < nrow_s1
                set(ax,'XTickLabel',[])
            end
            if shank > 1
                set(ax,'YTickLabel',[])
            end
        end
    end
    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = 'shape r (Fisher-z mean over trials)';
    xlabel(tl,'ms from stim onset (window centre)')
    ylabel(tl,'Hz')
    % title(tl,{[animal ' S1 probe, ' grpnames{2} ' (' num2str(ntr2) ' response trials): r across trials between V1 spike rhythmicity (responsive-channel mean) and each S1 channel''s power, per frequency and time'], ...
    %     ['red = trials with more S1 power there have more rhythmic V1 spiking there; ' num2str(tf_win_ms) ' ms windows; bold blue = CH' num2str(s1_refch)]},'Interpreter','none') %commented out with the across-trial map
    title(tl,{[animal ' S1 probe, ' grpnames{2} ' (' num2str(ntr2) ' response trials): shape match between V1 spike rhythmicity (responsive-channel mean) and each S1 channel''s power over time, per frequency'], ...
        ['r over a ' num2str(shape_block_ms) ' ms block of ' num2str(tf_win_ms) ' ms windows around each time, per trial, Fisher-z averaged; red = they rise and fall together; bold blue = CH' num2str(s1_refch)], ...
        ['outlined / not faded = cluster beats ' num2str(n_shuf) ' trial shuffles (pixel p < ' num2str(shuf_p) ', cluster p < ' num2str(clust_p) ', corrected over the whole S1 probe)']},'Interpreter','none')

    if plotsave == 1
        fname = sprintf([animal '-' 'S1_vs_V1_timefreq_rmaps_S1probe']);
        savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
    end
end
