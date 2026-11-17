% NNx_StochRes_MUA_singleanimal.m
%
% Description: first look for signs of stochastic resonance (SR) in V1 multi-unit activity (MUA), one
%   animal at a time. Compares two condition groups (e.g. light alone vs light + whisker, where the extra
%   input could act as the "noise"), each made of one or more conditions (subgroups, e.g. intensities):
%     1) MUA (0/1 spike onsets, 30 kHz) -> 1 ms spike counts -> Gaussian-smoothed spike rate per trial
%     1b) response trials: per condition, each trial is matched against a template of that condition's
%        typical response (leave-one-out trimmed mean, resp_ms, rebuilt from the passing trials tmpl_iter times);
%        a response when its match beats resp_null_prct of the same match at fake stim times in the pre-stim
%        baseline (mean z over the condition's responsive channels). Only
%        response trials are used from 2) on; the fraction of response trials per group is reported
%        (fewer misses with the added input would itself be an SR sign)
%     2) look at the trials: overlaid in the probe layout (how well they line up as is)
%     3) spike-rate spectrum (PSD) of each trial on the raw 1 ms counts, divided by the trial's mean rate
%        (a random/Poisson spike train is flat at 1, so rhythm = bumps above 1), condition means in the
%        probe layout
%     4) line the trials up by cross-correlation (one lag per trial, shared by all channels, capped at
%        +/- max_lag_ms) to the most typical response trial. Same alignment on the pre-stim baseline as a null
%        (alignment always raises correlation a bit, even on noise)
%     5) trial-to-trial correlation (as in NNx_LFPCorrelation_LvsLW.m), unaligned vs aligned:
%        does one group have more consistent trials once onset jitter is taken out?
%     6) S1 (probe 2) LFP spectra next to the V1 spike-rate spectra, when you say S1 is usable (do_S1):
%        if the whisker input acts as SR "noise", V1 spiking should get more rhythmic in the bands where
%        S1 gets more power, at the same times. For group 2 (e.g. LW), sliding-window spectra of V1 spike
%        rhythmicity and of every S1 channel's power per trial; at each frequency the shapes of the two time
%        courses are correlated over a sliding block, per trial, averaged, tested against trial shuffles
%        with probe-wide cluster correction (S1-probe layout of time-frequency maps)
%        (spike-LFP / LFP-LFP phase locking was tried and dropped, see commit 1747056)
%     7) post-run report (printed once at the end and saved) and the results saved with all settings
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
% Outputs: saved to <save_directory>\StochRes\<animal>\<animal>-StochRes_MUA.mat (-v7.3) in the last section,
%   results only (NOT spkcount, rate, lfp_s1: re-run the loading sections for those), plus:
%   settings (struct) - every Setup value of the run (+ load_ms, run_date), for pooling animals later
%   report (string, lines x 1) - the post-run summary, also printed at the end
% Workspace variables (and what's in the .mat):
%   prb (double) - probe used; pmap (double, rows x shanks) - its ProbeMaps raw channel IDs
%   chanlist (double, channels x 1) - kept raw channel IDs (probe minus Ch_Remove), sorted
%   chrow, chshank (double, channels x 1) - ProbeMaps row (1 = top) and shank of each kept channel
%   condis (cell, 1 x conditions) - condition names, group 1's first; condgrp (1 x conditions) - group 1 or 2
%   npair (double), pairs (double, pairs x 2) - condition pairs by position: [place in condis of the group 1
%     condition, of the group 2 condition]; pairnames (cell, 1 x pairs) - e.g. 'L_8 vs LW_8'. A group's
%     conditions are never pooled: groups are compared inside each pair
%   t_bin (double, 1 x bins) - ms (vs stim onset) where each count bin starts
%   spkcount (cell, 1 x conditions) - (single, trials x bins x channels) spike counts per bin_ms bin,
%     good trials only, files of a condition pooled, channels in chanlist order
%   trinfo (cell, 1 x conditions) - (double, trials x 2) file number within the condition, trial (tr_keep value)
%   ntrials (double, 1 x conditions) - good trials per condition
%   rate (cell, 1 x conditions) - (single, trials x bins x channels) spikes/s, spkcount smoothed with a
%     Gaussian kernel smooth_ms wide (total width)
%   filesused (string, files x 4) - condition, file name, good trials, trials left out by tr_conditional_use
%   (Quality Control: chan_medz_cond, respchans_cond, chan_medz, respscore, respnull, respthr, resp_tmpl and
%   tmpl_rounds are always made, but only saved with keep_qc = 1; resp, respchans, respfrac are always saved)
%   chan_medz_cond (single, channels x conditions) - each channel's median trial z: template-match score
%     (sum(trial .* template)/sum(template.^2) in resp_ms, trial minus its pre-stim median, template = trimmed
%     leave-one-out mean of the condition's other (response) trials) vs that score at the fake no-stim times
%   respchans_cond (logical, channels x conditions) - true = responsive in that condition (z >= chan_z_min)
%   chan_medz (single, channels x 1) - best chan_medz_cond over the conditions; respchans (logical,
%     channels x 1) - responsive in at least one condition (used by the later sections)
%   respscore (cell, 1 x conditions) - (single, trials x 1) mean z over the condition's responsive channels
%   resp (cell, 1 x conditions) - (logical, trials x 1) true = response trial (respscore > respthr).
%     Only these trials are used from here on
%   respnull (cell, 1 x conditions) - (single, trials*fake windows x 1) the same score at the fake no-stim
%     times; respthr (double, 1 x conditions) - its resp_null_prct percentile (NaN = no responsive channel)
%   resp_tmpl (cell, 1 x conditions) - (single, resp_ms bins x channels) spikes/s above pre-stim, final
%     template (trimmed mean of the response trials); tmpl_rounds (double, 1 x conditions) - rounds run
%   respfrac (double, 1 x conditions) - fraction of response trials
%   rate_mean, rate_sem (cell, 1 x conditions) - (single, bins x channels) spikes/s, mean and SEM of
%     rate over the response trials
%   psd_f (double, freqs x 1) - Hz of the spectra (psd_hz, multitaper nfft steps)
%   psd_tot (cell, 1 x conditions) - (single, response trials x freqs x channels) multitaper spectrum of
%     the raw counts in psd_ms divided by the trial's random (Poisson) firing level 2*mean count/fs:
%     1 = no rhythm. NaN at line_hz if line_out = 1
%   psdm_tot (cell, 1 x conditions) - (double, freqs x channels) mean over trials
%   (do_align = 1) align_methods (cell, 1 x 1) - 'typical' (template = the response trial whose rate, summed
%     over the responsive channels, correlates best with the mean of the other response trials in align_use_ms)
%   align_ref (double, 1 x conditions) - that reference trial (row of trinfo{con}); align_ref_r (double, 1 x
%     conditions) - its r with the others' mean. Only saved / reported with keep_qc = 1
%   lags, lags_null (cell, 1 x methods){1 x conditions} - (double, response trials x 1) ms, one shift per
%     trial from cross-correlating the responsive-channel summed rate in align_ms (null: same-length
%     pre-stim window) with the template, |lag| <= max_lag_ms. Aligned trial(t) = trial(t + lag)
%   align_use_ms (double, 1 x 2) - ms, the alignment window actually used: align_ms trimmed by max_lag_ms
%     where it would run off the loaded data (load_ms); the null window is the same length, ending at 0
%   lagpk (cell, 1 x methods){1 x conditions} - (double, response trials x 1) peak cross-correlation r
%   lagsd, lagsd_null (double, methods x conditions) - ms, SD of the lags (trial-to-trial jitter)
%   (do_align = 1, aligned-trial section) al_names (cell, 1 x versions) - 'unaligned', 'typical'
%   tc_trial (cell, 1 x versions){1 x conditions} - (double, response trials x channels) each trial's mean r
%     (smoothed rate time courses in align_use_ms, Fisher z mean, back to r) with the other trials of its condition
%   tc_cond, tc_cond_null (cell, 1 x versions) - (double, channels x conditions) mean pairwise trial r within
%     each condition; null = the same on the pre-stim window aligned with the null lags
%   tc_p (cell, 1 x versions) - (double, channels x pairs) label-shuffle permutation p, group 2 vs group 1
%     condition of each pair (only that pair's trials shuffled)
%   (do_S1 = 1) lfp_s1 (cell, 1 x conditions) - (single, response trials x samples x S1 channels) uV, S1 LFP
%     in lfp_load_ms (1 kHz), same trials/order as the MUA; t_lfp (double, 1 x samples) - ms of each sample
%   chanlist_s1 (double, channels x 1) - kept S1 raw channel IDs; prb_s1 - S1 probe number
%   s1_refch (double) - raw channel ID of the S1 reference channel (kept channel s1_ref_row of the peak
%     s1_shank_label shank)
%   (do_S1 = 1, time-frequency, group 2's response trials only, EACH group 2 condition on its own: the cells
%   below are {1 x conditions}, empty for group 1 and for group 2 conditions with < 3 response trials)
%   tf_t (double, 1 x windows) - ms, sliding-window centres used (tf_ms, trimmed to the loaded data);
%     tf_f (double, freqs x 1) - Hz (psd_hz, steps of 1000/nfft of tf_win_ms)
%   v1tr {con} (double, freqs x windows x trials) - V1 spike rhythmicity per response trial: multitaper
%     power of the raw counts per window / that window's random-firing level, responsive-channel mean, dB
%   ntr_s1 (double, 1 x conditions) - response trials used per condition (0 = not analysed)
%   r_shape_ch {con} (double, freqs x windows x S1 channels) - shape match: per trial, r over time (a shape_block_ms
%     block of windows centred on each window) between V1 rhythmicity and that S1 channel's power at each
%     frequency, Fisher-z averaged over the condition's trials (NaN where the block runs off either end). Plotted
%   shape_zpix, shape_ppix {con} (double, freqs x windows x S1 channels) - each pixel's z and two-sided p against
%     n_shuf trial shuffles within the condition (V1 of one trial paired with another trial's S1)
%   shape_sig {con} (logical, freqs x windows x S1 channels) - pixels in clusters (|z| past shuf_p, positive and
%     negative apart) whose summed |z| beats shape_clustthr(con) = the (1 - clust_p) percentile of the biggest
%     chance cluster per shuffle over all S1 channels (per condition)
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
%one condition per group is fine. More than one (e.g. intensities): give both groups the same number, IN
%MATCHING ORDER, because they are paired by position (1st with 1st: L_4 vs LW_4, 2nd with 2nd: L_8 vs LW_8).
%Conditions are never pooled within a group; the groups are compared inside each pair
grpnames = {'L','LW'};
condis_grp{1} = {'L_4','L_8','L_12','L_15'}; %conditions in group 1
condis_grp{2} = {'LW_4','LW_8', 'LW_12', 'LW_15'}; %conditions in group 2

%leave out trials flagged under these conditional trial set names? (made with NNx_TrialRemove_Conditional.m,
%e.g. {"whisker_twitch_artifact"}). {} = no, only tr_remove
tr_conditional_use = {};

%where are the spiking, LFP and ProbeInfo files?
spike_dir = 'E:\Roy\Processed Silicon Probe Data\Spiking';
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%Where should I, Matlab, put the results? (Desktop while we're testing)
%save_directory = 'C:\Users\rmissa\Desktop\StochRes_test'; %results in \StochRes\<animal>, figures in \Figures_StochRes\<animal>
save_directory = 'E:\Roy\Processed Silicon Probe Data\StochRes';

%%%Analysis details%%%
stim_onset_30k = 150001; %sample of stimulus onset (0 ms) in stim_spike_stimchunks, 30 kHz
stim_onset_lfp = 5001; %sample of stimulus onset (0 ms) in stim_lfp_stimchunks, 1 kHz
fs_mua = 30000; %Hz, MUA sampling rate
win_ms = [0 3000]; %ms post stimulus, the response window everything is measured in
% base_ms = [-3000 0]; %commented out: the 2000 ms response template needs more pre-stim room for its fake windows
base_ms = [-5000 0]; %ms, pre-stim stretch loaded (fake windows for the response trials + alignment null). The files hold 5000 ms pre-stim, no more
bin_ms = 1; %ms, spike count bin size (1 ms -> spike rate and spectra at 1 kHz)
smooth_ms = 25; %ms, total width of the Gaussian kernel for the spike-rate traces (not the SD; ~ SD x 5)

%%%Response trial details%%%
%only trials with a response are used from the trial plots on (no-response trials just muddy things).
%template matching, per condition: a trial's score = how much of its condition's typical response shape
%(the template: trimmed mean of the OTHER trials, resp_ms, each minus its pre-stim median) it holds, as a z
%against the same score at fake stim times in the pre-stim baseline (fake_ms), averaged over that
%condition's responsive channels. The template is then rebuilt from the trials that passed (tmpl_iter rounds)
%(the old rate-rise score, 0-300 ms minus -200-0 ms, missed weak low-light responses and caught bursts.
%Before that, tried and dropped: resp_ms [75 500], [500 1000], [500 1500], [75 300] with the "before"
%stretch ending at 70 ms; pre_resp_ms 500; median/MAD z)
% resp_ms = [0 300]; %commented out: rate-rise window, the template covers the whole response now
% resp_ms = [0 1000]; %commented out: widened to hold the whole response
resp_ms = [0 2000]; %ms post stimulus, the stretch the template covers (the full response, P1 + P2)
pre_resp_ms = 200; %ms, length of the "before" stretch whose median is each trial's baseline
pre_end_ms = 0; %ms vs the (real or fake) stimulus, where the "before" stretch ends (0 -> -200 to 0 ms, just before the stimulus)
% fake_ms = -2700:300:-300; %commented out: fit the 300 ms window, the 1000 ms window needs earlier fakes
% fake_ms = -2800:150:-1000; %commented out: fit the 1000 ms resp_ms
% fake_ms = -2800:100:-2000; %commented out: all that fit with 3000 ms pre-stim loaded, and each shared 1900 of its 2000 ms with the next
fake_ms = -4800:250:-2000; %ms, fake stimulus times in the pre-stim baseline (no-stim reference, 12 per trial, neighbours share 1750 ms). Needs fake_ms(1) + pre_end_ms - pre_resp_ms >= base_ms(1) and fake_ms(end) + resp_ms(2) <= 0
tmpl_trim = 20; %percent of trials left out (half each end) at every time point when averaging the template, so one trial's burst can't make a bump in it
tmpl_iter = 2; %rounds of rebuilding the template from the trials that passed (stops earlier if the calls stop changing)
chan_z_min = 1; %a channel counts as responsive (in a condition) when its typical (median) trial z is at least this
resp_null_prct = 99; %a trial is a response trial when its score beats this percentile of its condition's fake scores

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
    % 1 = yes, to the most typical response trial (best shape match with the others in align_ms)
align_ms = [350 3000]; %ms post stimulus, the stretch the cross-correlation looks at (sustained/rhythmic part)
max_lag_ms = 50; %ms, biggest shift allowed for a trial (either way)
% woody_iter = 3; %commented out: the leave-one-out (Woody) template was dropped, a single reference trial is used
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

%%%Quality Control%%%
%want to dig into how the early steps handled the data (which trials were kept, which channels respond, what
%the response detection compared against)? Everything for that is switched on/off here
%responsive channels on the probe? (each channel's best typical z over the conditions, responsive ones circled)
plot_respchans = 0;
    % 0 = no
    % 1 = yes (one probe-map figure, end of the Response trials section)
%every trial of each condition as a heatmap, with which got called a response + mean response vs no-response traces,
%plus each condition's response template?
plot_respcheck = 0;
    % 0 = no
    % 1 = yes (one figure per condition + one template figure)
%save the Quality Control figures? (own subfolder: <save_directory>\Figures_StochRes\<animal>\QC, apart from the
%results figures; separate from plotsave)
qc_plotsave = 0;
    % 0 = no, just look at them
    % 1 = yes
%keep the response detection details? (saved with the results + listed in the post-run report)
keep_qc = 0;
    % 0 = no: only the trial calls (resp), responsive channels (respchans) and response fractions are kept
    % 1 = yes: also per-trial scores, the no-stim scores + thresholds, per-condition channel z and responsive
    %     channels, the templates and template rounds, and the alignment reference trial of each condition (+ its
    %     shape r); the report lists every condition's kept (response) and left-out trials, its responsive
    %     channels and its alignment reference trial

%%%Plotting details%%%
%overlaid single trials?
plot_trials = 1;
    % 0 = no
    % 1 = yes (probe layout, one figure per condition)
% stack_offset = 100; %commented out: the stacked-trial plots (scroll viewer, every-trial check) were dropped
plot_ms = [-500 3000]; %ms, time range shown in the trial plots
%probe-layout spectra?
plot_psd = 1;
%alignment method comparison?
plot_align = 1;
%trial-to-trial correlation maps?
plot_corr = 1;
%S1 LFP spectra vs V1 spike-rate spectra? (do_S1 = 1)
plot_lfpcoh = 1;
%save the figures? (as PNGs; a figure with a slider gets one PNG per panel, e.g. ...-psd_tot-L_8_vs_LW_8.png)
plotsave = 1;
png_res = 150; %dpi of the saved PNGs (a full-screen figure -> ~2900 x 1600 px)

%delete big in-between variables as we go? (saves memory)
mr_clean = 1;

%%%Speed%%%
%run the heaviest loop (the S1 shape-match shuffles, one S1 channel per worker) on several CPU cores?
%Needs the Parallel Computing Toolbox; without it everything runs on one core as normal (MATLAB's own maths
%still uses all cores). Starting the worker pool takes ~20-60 s the first time in a MATLAB session
n_workers = Inf;
    % 0 = one core
    % Inf = as many workers as the pool has (or a number, e.g. 4)

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

%condition pairs: the k-th condition of group 1 goes with the k-th of group 2 (e.g. L_8 with LW_8).
%Conditions are never pooled within a group; the groups are compared inside each pair
if numel(condis_grp) ~= 2 || numel(condis_grp{1}) ~= numel(condis_grp{2})
    error('Conditions are paired by position (condis_grp{1}{k} with condis_grp{2}{k}): give 2 groups with the same number of conditions')
end
npair = numel(condis_grp{1}); %number of condition pairs
pairs = [(1:npair)' numel(condis_grp{1}) + (1:npair)']; % [pairs x 2] place in condis of the group 1 and group 2 condition
pairnames = cell(1,npair); %{1 x pairs} e.g. 'L_8 vs LW_8'
for k = 1:npair
    pairnames{k} = [condis{pairs(k,1)} ' vs ' condis{pairs(k,2)}];
end

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

%Did the trial respond at all? Template matching, one condition at a time (each condition gets its own
%picture of a typical response, so dim light is compared with dim light, not with bright light):
%  1) every trial's smoothed rate in resp_ms, minus that trial's median rate in the "before" stretch
%     (a median, so a spontaneous burst there barely moves it)
%  2) each trial's template = trimmed mean of the OTHER trials (leave-one-out, so a trial can't match
%     itself; trimmed, so one trial's burst can't make a bump in it). Round 1 uses all the trials
%  3) score = sum(trial .* template) / sum(template.^2): how many templates' worth of that shape the
%     trial holds (1 = a full one, 0 = nothing template-shaped). The same score at fake stimulus times in
%     the pre-stim baseline (fake_ms, against the same template) says what "nothing happened" looks like
%  4) per channel, scores in units of the no-stim spread (z); responsive channels (this condition) =
%     median trial z >= chan_z_min; trial score = mean z over them; response = beats resp_null_prct of
%     the fake scores
%  5) rebuild the templates from the trials that passed (sharper picture) and redo 3-4, tmpl_iter times
%     at most, stopping early once the calls stop changing
%Only the shape decides: a faint template (few responders) shrinks the real and the fake scores alike.
%Uses the smoothed rate (template from 1 ms counts would be mostly noise); the 25 ms kernel leaks only a
%few ms of response into the "before" stretch, and its median ignores that
disp('MUA: Finding the trials with a response (template matching)')
tic

%do the fake windows fit? (inside the loaded stretch, and ending before the real stimulus)
if fake_ms(1) + pre_end_ms - pre_resp_ms < load_ms(1)
    error(['The first fake window starts before the loaded data (' num2str(load_ms(1)) ' ms). Move fake_ms later or base_ms earlier'])
end
if fake_ms(end) + resp_ms(2) > 0
    error('The last fake response window runs past the real stimulus. Move fake_ms earlier')
end

respidx = find(t_bin >= resp_ms(1) & t_bin < resp_ms(2)); %bins the template covers
preidx = find(t_bin >= pre_end_ms - pre_resp_ms & t_bin < pre_end_ms); %bins of the "before" stretch
nw = numel(respidx); %template length in bins
nfake = numel(fake_ms); %fake windows per trial

respscore = cell(1,ncond); %{conditions} [trials x 1] trial score (mean z over the condition's responsive channels)
resp = cell(1,ncond); %{conditions} [trials x 1] true = response trial
respnull = cell(1,ncond); %{conditions} [(trials*fake windows) x 1] the same score at the fake no-stim times
respthr = nan(1,ncond); %its resp_null_prct percentile, per condition
chan_medz_cond = nan(nch,ncond); % [channels x conditions] each channel's median trial z
respchans_cond = false(nch,ncond); % [channels x conditions] true = responsive in that condition
resp_tmpl = cell(1,ncond); %{conditions} [template bins x channels] spikes/s above baseline, final template
tmpl_rounds = zeros(1,ncond); %template rounds actually run per condition
for con = 1:ncond
    ntr = ntrials(con);
    disp(['MUA: Template matching ' condis{con}])

    %real windows, each trial minus its own pre-stim median
    X = rate{con}(:,respidx,:) - median(rate{con}(:,preidx,:),2); % [trials x template bins x channels] spikes/s

    %fake windows, the same way
    Xf = zeros(ntr,nw,nch,nfake,'single'); % [trials x template bins x channels x fake windows]
    for f = 1:nfake
        fresp = find(t_bin >= fake_ms(f)+resp_ms(1) & t_bin < fake_ms(f)+resp_ms(2)); %fake response window
        fpre = find(t_bin >= fake_ms(f)+pre_end_ms-pre_resp_ms & t_bin < fake_ms(f)+pre_end_ms); %and its "before" stretch
        Xf(:,:,:,f) = rate{con}(:,fresp,:) - median(rate{con}(:,fpre,:),2);
    end

    use = true(ntr,1); %trials the templates are built from (round 1: all of them)
    prevresp = []; %last round's calls
    for rnd = 1:tmpl_iter+1
        %each trial's template: trimmed mean of the other trials in use
        T = zeros(ntr,nw,nch,'single'); % [trials x template bins x channels]
        for tr = 1:ntr
            others = use;
            others(tr) = false;
            T(tr,:,:) = trimmean(X(others,:,:), tmpl_trim, 1);
        end
        den = sum(T.^2,2); % [trials x 1 x channels] size of each template
        den(den == 0) = NaN; %flat template (silent channel): no score

        %scores, real and fake
        sc_ch = reshape(sum(X.*T,2)./den, ntr, nch); % [trials x channels]
        scf_ch = zeros(ntr,nfake,nch,'single');
        for f = 1:nfake
            scf_ch(:,f,:) = sum(Xf(:,:,:,f).*T,2)./den; %each fake window against the same trial's template
        end
        scf_ch = reshape(scf_ch, ntr*nfake, nch); % [(trials*fake windows) x channels]

        %each channel's scores in units of its no-stim spread
        fakemu = mean(scf_ch,1,'omitnan'); % [1 x channels]
        fakesd = std(scf_ch,0,1,'omitnan'); % [1 x channels]
        fakesd(fakesd == 0) = NaN; %silent channel, drops out
        zch = (sc_ch - fakemu) ./ fakesd; % [trials x channels]
        zfake = (scf_ch - fakemu) ./ fakesd; % [(trials*fake windows) x channels]

        %this condition's responsive channels
        chan_medz_cond(:,con) = median(zch,1,'omitnan')';
        respchans_cond(:,con) = chan_medz_cond(:,con) >= chan_z_min;
        tmpl_rounds(con) = rnd;
        if any(respchans_cond(:,con)) == false
            respscore{con} = nan(ntr,1);
            resp{con} = false(ntr,1);
            respnull{con} = nan(ntr*nfake,1);
            respthr(con) = NaN;
            run_notes(end+1) = string(condis{con}) + ": no channel with a typical trial z >= chan_z_min in template round " + rnd + ", no response trials called";
            break
        end

        %trial scores and the threshold
        respscore{con} = mean(zch(:,respchans_cond(:,con)),2); % [trials x 1]
        respnull{con} = mean(zfake(:,respchans_cond(:,con)),2); % [(trials*fake windows) x 1]
        respthr(con) = prctile(respnull{con}, resp_null_prct); %a response has to beat this
        resp{con} = respscore{con} > respthr(con);

        %rebuild the templates from the trials that passed? (not if nothing changed or too few passed)
        if isequal(resp{con}, prevresp)
            break
        end
        if sum(resp{con}) < 3
            run_notes(end+1) = string(condis{con}) + ": only " + sum(resp{con}) + " response trial(s) in template round " + rnd + ", the template wasn't rebuilt from them";
            break
        end
        prevresp = resp{con};
        use = resp{con};
    end

    %the final template, for the pictures (trimmed mean of the response trials, no leave-one-out)
    if sum(resp{con}) >= 2
        resp_tmpl{con} = reshape(trimmean(X(resp{con},:,:), tmpl_trim, 1), nw, nch);
    else
        resp_tmpl{con} = nan(nw,nch,'single');
    end
end
clear X Xf T den sc_ch scf_ch zch zfake

%channels responsive in any condition (what the later sections use: alignment, plot labels)
chan_medz = max(chan_medz_cond,[],2); % [channels x 1] best typical trial z over the conditions
respchans = any(respchans_cond,2); % [channels x 1]
if any(respchans) == false
    error(['No channel has a typical trial z >= chan_z_min (' num2str(chan_z_min) ') in any condition, so no responses can be called. Try a lower chan_z_min'])
end

%{
%commented out: the old rate-rise score (0-300 ms minus -200-0 ms, raw counts, conditions pooled). Missed weak
%low-light responses and called spontaneous bursts responses. Template matching above instead
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
row = 0;
for con = 1:ncond
    idx = row + (1:ntrials(con));
    respscore{con} = sc(idx);
    resp{con} = sc(idx) > respthr;
    row = row + ntrials(con);
end
clear allrise allfake zfake sc
%}

%fraction of response trials per condition
respfrac = zeros(1,ncond);
for con = 1:ncond
    respfrac(con) = mean(resp{con});
end

%commented out: pooled a group's conditions. Conditions stay apart, the report compares them per pair
% respfrac_grp = zeros(1,numel(condis_grp));
% for g = 1:numel(condis_grp)
%     respfrac_grp(g) = sum(cat(1, resp{condgrp == g})) / sum(ntrials(condgrp == g));
% end

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

elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);

%(the temporary check section that followed, with the score-vs-no-stim histograms and the template
%picture, was dropped once the template matching was good enough; the every-trial heatmaps and the
%template picture are in their own optional section below, plot_respcheck)

%responsive channels on the probe (Quality Control: plot_respchans): each channel's best typical z over the
%conditions, responsive ones circled
if plot_respchans == 1
    qc_dir = [save_directory '\Figures_StochRes\' animal '\QC']; %Quality Control figures go here, apart from the results figures
    if qc_plotsave == 1
        %create QC figure folder
        if isfolder(qc_dir) == false
            mkdir(qc_dir)
        end
    end

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
    title([animal ' ' ProbeInfo.Areas{prb} '-Probe: best typical response z over the conditions (o = responsive in at least one, z >= ' num2str(chan_z_min) '; x = removed)'])

    if qc_plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'resptrials_channels']);
        % savefig(gcf,[qc_dir '\' fname],'compact') %commented out: PNGs for now
        exportgraphics(gcf,[qc_dir '\' fname '.png'],'Resolution',png_res)
    end
end


%% Data Visualization: Every trial, response or not
%(Quality Control: plot_respcheck)

%per condition: every trial as a heatmap (rate summed over the condition's responsive channels, recording
%order), a strip beside it marking which got called a response, and the mean of the response vs
%no-response trials underneath. Then one figure of every condition's final response template
if plot_respcheck == 1
    disp('MUA: Every-trial heatmaps, by response call')

    qc_dir = [save_directory '\Figures_StochRes\' animal '\QC']; %Quality Control figures go here, apart from the results figures
    if qc_plotsave == 1
        %create QC figure folder
        if isfolder(qc_dir) == false
            mkdir(qc_dir)
        end
    end

    pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
    t_plot = t_bin(pidx); %ms

    for con = 1:ncond
        %trsum = sum(rate{con}(:,pidx,respchans),3); %commented out: responsive channels are per condition now
        trsum = sum(rate{con}(:,pidx,respchans_cond(:,con)),3); % [trials x bins shown] rate summed over this condition's responsive channels
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

        if qc_plotsave == 1
            fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-' 'resptrials_trials']);
            % savefig(gcf,[qc_dir '\' fname],'compact') %commented out: PNGs for now
            exportgraphics(gcf,[qc_dir '\' fname '.png'],'Resolution',png_res)
        end
    end

    %each condition's final template (what a typical response looks like), summed over its responsive channels
    grpcol = [0 0.45 0.74; 0.85 0.33 0.10]; %colour of group 1, group 2 (shaded darker -> lighter by condition)
    t_tmpl = t_bin(respidx); %ms
    figure()
    hold on
    lnames = {};
    for con = 1:ncond
        if any(respchans_cond(:,con)) == false || all(isnan(resp_tmpl{con}(:)))
            continue
        end
        plot(t_tmpl, sum(resp_tmpl{con}(:,respchans_cond(:,con)),2), 'LineWidth',1.5, 'Color',grpcol(condgrp(con),:)*(0.5 + 0.5*con/ncond))
        lnames{end+1} = [condis{con} ' (' num2str(sum(respchans_cond(:,con))) ' ch, ' num2str(sum(resp{con})) ' trials)'];
    end
    xline(0)
    yline(0,':')
    hold off
    xlabel('ms from stim onset')
    ylabel('spikes/s above pre-stim, summed over responsive channels')
    legend(lnames,'Location','northeast','Interpreter','none')
    title([animal ' ' ProbeInfo.Areas{prb} '-Probe: final response template per condition'])

    if qc_plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'resptrials_templates']);
        % savefig(gcf,[qc_dir '\' fname],'compact') %commented out: PNGs for now
        exportgraphics(gcf,[qc_dir '\' fname '.png'],'Resolution',png_res)
    end
end



%Now beginning the single trial pictures (response trials only from here on)

%% Data Visualization: Overlaid trials

%the response trials all on top of each other with nothing shifted (grey), plus the condition mean
%(black), every channel in the probe layout. How well do the trials line up as they are?
%ONE figure for all conditions: each condition is drawn in its own panel, and the slider at the bottom
%shows one panel at a time (click / drag it, or click it and use the arrow keys). Everything is drawn
%up front, so the script doesn't wait here: it carries on while you scroll.
%(was one figure per condition)
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

    ovfig = figure();
    set(ovfig, 'Position', get(0, 'Screensize'));
    for con = 1:ncond
        %this condition's panel (only the first one starts visible); Tag + UserData let the slider find it
        pan = uipanel(ovfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
            'Tag','condpanel', 'UserData',con, 'Visible', con == 1);

        rtrs = find(resp{con} == true); %response trials
        nresp = numel(rtrs);
        if nresp == 0
            run_notes(end+1) = string(condis{con}) + ": no response trials, no overlaid trial picture";
            ax = axes(pan);
            axis(ax,'off')
            text(ax,0.5,0.5,[condis{con} ': no response trials'],'Units','normalized','HorizontalAlignment','center','Interpreter','none')
            continue
        end
        trdat = rate{con}(rtrs,pidx,:); % [response trials x bins shown x channels] spikes/s
        ytop = max(trdat,[],'all') * 1.05; %one y range for the whole probe

        % figure() %commented out: one figure per condition before, now one panel per condition
        % set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout(pan,nrow,nshank,'TileSpacing','none','Padding','tight');
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
        title(tl,[animal ' ' ProbeInfo.Areas{prb} '-Probe ' condis{con} ' (condition ' num2str(con) ' of ' num2str(ncond) '): overlaid response trials (grey, ' num2str(nresp) ' of ' num2str(ntrials(con)) ') and their mean (black)'],'Interpreter','none')

        % if plotsave == 1 %commented out: one file per condition before, now one file for all (below)
        %     fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-' 'trials_overlaid']);
        %     savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact')
        % end
    end

    %the slider: moving it shows the panel whose UserData (condition number) matches the slider value and hides
    %the rest. It finds the panels by their Tag, so it still works when the saved figure is opened again later
    if ncond > 1
        uicontrol(ovfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
            'String','condition:', 'HorizontalAlignment','right')
        uicontrol(ovfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
            'Min',1, 'Max',ncond, 'Value',1, 'SliderStep',[1 1]/(ncond-1), ...
            'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
            num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
    end
    drawnow %show it now; the script doesn't wait, it carries on

    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'trials_overlaid']);
        % savefig(ovfig,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now, one per panel
        %one PNG per panel (= condition): show only panel k, save its layout (not the slider), then back to panel 1
        ptags = condis; %file name tag of each panel
        pans = findobj(ovfig,'Tag','condpanel');
        for k = 1:numel(ptags)
            set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
            pk = findobj(pans,'flat','UserData',k);
            tlk = findobj(pk,'Type','tiledlayout');
            if isempty(tlk)
                tlk = pk; %no layout (just a note): the panel itself
            end
            exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' ptags{k} '.png'],'Resolution',png_res)
        end
        set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
        set(findobj(ovfig,'Style','slider'),'Value',1);
    end
end




%Now beginning the spectrum section

%% Data Analysis: Spike-rate spectrum (PSD) per trial and condition mean

%Is the firing rhythmic, and at what frequency? Multitaper spectrum (pmtm) of each response trial's raw
%1 ms spike counts in psd_ms (the sustained response, after the onset burst), per channel.
%Each trial's spectrum is divided by what random (Poisson) firing at that trial's own rate would give,
%so 1 = no rhythm, > 1 = rhythm at that frequency, and conditions with different rates compare fairly.
%(the pre-stim, evoked and induced versions were tried and dropped)
disp('MUA: Spike-rate spectra')
tic

fs = 1000/bin_ms; %Hz, sampling rate of the counts
psdidx = find(t_bin >= psd_ms(1) & t_bin < psd_ms(2)); %bins of the response window
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

psd_tot = cell(1,ncond); %{conditions} [response trials x freqs x channels] normalized (1 = random firing)
psdm_tot = cell(1,ncond); %{conditions} [freqs x channels] mean over trials
for con = 1:ncond
    rtrs = find(resp{con} == true); %response trials
    nresp = numel(rtrs);
    psd_tot{con} = nan(nresp,nf,nch,'single');
    if nresp < 2
        run_notes(end+1) = string(condis{con}) + ": fewer than 2 response trials, no spectra";
        psdm_tot{con} = nan(nf,nch);
        continue
    end

    for ch = 1:nch
        X = double(reshape(spkcount{con}(rtrs,psdidx,ch), nresp, nwin))'; % [samples x trials] counts, response window
        m = mean(X,1); % [1 x trials] mean count per bin of each trial
        poiss = 2*m/fs; %spectrum level of random (Poisson) firing at each trial's rate (one-sided)
        poiss(poiss == 0) = NaN; %a trial without spikes has no spectrum

        %each trial minus its own mean
        ptot = pmtm(X - m, psd_nw, nfft, fs); % [nfft/2+1 x trials], one column per trial
        psd_tot{con}(:,:,ch) = (ptot(fkeep,:) ./ poiss)';
    end

    psd_tot{con}(:,linemask,:) = NaN; %line noise out (if asked)
    psdm_tot{con} = reshape(mean(psd_tot{con},1,'omitnan'), nf, nch); %condition mean of the trial spectra
end
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data Visualization: Spike-rate spectra in the probe layout

%the spike-rate spectrum on every channel in the probe layout, one condition pair at a time (slider; the
%last panel has every condition). Dashed line at 1 = random firing, bumps above it = rhythm
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

    %ONE figure with a slider: a panel per condition pair (its two conditions, group 1 dark blue, group 2
    %dark red), then a last panel with every condition (blues / reds, dark to light). Everything is drawn up
    %front, so the script doesn't wait here. (was one figure with every condition in each tile; the
    %pre-stim / evoked / induced spectra were dropped before that)
    panelcons = cell(1,npair+1); %{panels} conditions shown in each panel
    paneltitle = cell(1,npair+1);
    for k = 1:npair
        panelcons{k} = pairs(k,:);
        paneltitle{k} = [pairnames{k} ' (pair ' num2str(k) ' of ' num2str(npair) ')'];
    end
    panelcons{npair+1} = 1:ncond;
    paneltitle{npair+1} = 'all conditions';
    npanel = npair + 1;

    %one y range for every panel, so they compare
    ytop = 0;
    for con = 1:ncond
        ytop = max(ytop, max(psdm_tot{con},[],'all','omitnan'));
    end
    if isnan(ytop) || ytop == 0
        ytop = 1;
    end
    ytop = ytop*1.05;

    psdfig = figure();
    set(psdfig, 'Position', get(0, 'Screensize'));
    for pnl = 1:npanel
        %this panel (only the first starts visible); Tag + UserData let the slider find it
        pan = uipanel(psdfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
            'Tag','condpanel', 'UserData',pnl, 'Visible', pnl == 1);
        cset = panelcons{pnl}; %its conditions
        if pnl <= npair
            pcol = grpdark(condgrp(cset),:); %a pair: one dark colour per group
        else
            pcol = condcol(cset,:); %everything: group shades
        end

        tl = tiledlayout(pan,nrow,nshank,'TileSpacing','none','Padding','tight');
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
                for i = 1:numel(cset)
                    plot(ax, psd_f, psdm_tot{cset(i)}(:,ch), 'Color',pcol(i,:), 'LineWidth',1)
                end
                yline(ax,1,'k--') %random firing = 1
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
        %legend in its own strip under the probe, never on top of the data (only the condition lines, not the dashed line)
        lns = findobj(axlast,'Type','line','LineStyle','-');
        lgd = legend(axlast, flipud(lns), legnames(cset), 'Orientation','horizontal', 'Interpreter','none');
        lgd.Layout.Tile = 'south';
        xlabel(tl,'Hz')
        ylabel(tl,'power / random firing (1 = no rhythm)')
        title(tl,{[animal ' ' ProbeInfo.Areas{prb} '-Probe spike-rate spectrum (' num2str(psd_ms(1)) '-' num2str(psd_ms(2)) ' ms): ' paneltitle{pnl}], ...
            ['blue = ' grpnames{1} ', red = ' grpnames{2} ', dashed line = random firing (1), bold blue CH = responsive channel']},'Interpreter','none')
    end

    %the slider (finds the panels by their Tag, so it still works in the saved figure)
    uicontrol(psdfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
        'String','pair (last = all):', 'HorizontalAlignment','right')
    uicontrol(psdfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
        'Min',1, 'Max',npanel, 'Value',1, 'SliderStep',[1 1]/(npanel-1), ...
        'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
        num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
    drawnow %show it now; the script doesn't wait, it carries on

    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'psd_tot']);
        % savefig(psdfig,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now, one per panel
        %one PNG per panel (= pair, then all conditions): show only panel k, save its layout (not the slider), then back to panel 1
        ptags = [strrep(pairnames,' ','_') {'allconds'}]; %file name tag of each panel, e.g. L_8_vs_LW_8
        pans = findobj(psdfig,'Tag','condpanel');
        for k = 1:numel(ptags)
            set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
            pk = findobj(pans,'flat','UserData',k);
            tlk = findobj(pk,'Type','tiledlayout');
            exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' ptags{k} '.png'],'Resolution',png_res)
        end
        set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
        set(findobj(psdfig,'Style','slider'),'Value',1);
    end
end




%Now beginning the trial alignment section

%% Data Analysis: Trial alignment by cross-correlation

%Maybe the trials share a rhythm but start a little earlier or later. Each response trial gets ONE shift
%(lag, the same for every channel), found by cross-correlating its spike rate (summed over the responsive
%channels, in align_ms) with a template, biggest shift max_lag_ms either way. Template = the most typical
%response trial ('typical'): the one whose rate in the alignment window correlates best with the mean of the
%other response trials, so it's picked by shape over the stretch being aligned, not by how big it is.
%(the highest-respscore trial was used before: that score grows with response size and covers resp_ms, not
%align_ms. The leave-one-out and first-trial templates were tried and dropped before that)
%Null: the exact same thing on a same-length stretch before the stimulus, where nothing should line up.
%Any alignment makes trials look a bit more alike, even pure noise, so the later sections compare to this.
%Lag convention: a lag of +d ms means the trial ran d ms late, so the aligned trial(t) = trial(t + d)
if do_align == 1
    disp('MUA: Lining the trials up')
    tic

    % align_methods = {'loo','strongest','first'}; %commented out: only the strongest-trial template is kept
    % align_methods = {'strongest'}; %commented out: reference now picked by shape, not by the response score
    align_methods = {'typical'}; %(kept as a list so the results keep their {method}{condition} layout)
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

    lags = cell(1,nmeth); %{methods}{conditions} [response trials x 1] ms, real window
    lags_null = cell(1,nmeth); %{methods}{conditions} [response trials x 1] ms, pre-stim window
    lagpk = cell(1,nmeth); %{methods}{conditions} [response trials x 1] peak cross-correlation (r) at that lag
    lagsd = nan(nmeth,ncond); %ms, SD of the lags (jitter), real window
    lagsd_null = nan(nmeth,ncond); %ms, the same on the pre-stim window
    align_ref = nan(1,ncond); %reference trial of each condition (row of trinfo{con})
    align_ref_r = nan(1,ncond); %its shape match: r with the mean of the other response trials in the alignment window
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
        % [~,strongest] = max(respscore{con}(rtrs)); %commented out: the highest score favours big responses, see above

        %reference trial: best shape match (r) with the mean of the OTHER response trials, in the real alignment window
        Xw = sig(:,winbins{1}); % [response trials x window bins]
        shaper = nan(nresp,1); % [response trials x 1] r of each trial with the others' mean
        for tr = 1:nresp
            othermean = mean(Xw([1:tr-1 tr+1:nresp],:),1); % [1 x window bins]
            shaper(tr) = corr(Xw(tr,:)', othermean'); %NaN for a flat trial (no spikes), max skips it
        end
        [align_ref_r(con),reftr] = max(shaper); %the most typical trial
        align_ref(con) = rtrs(reftr); %its row in trinfo{con}

        for w = 1:2 %1 = real window, 2 = null (pre-stim) window (same reference trial for both)
            wb = winbins{w}; %bins of this window
            trlag = zeros(nresp,1); %bins, lag of each trial
            trpk = nan(nresp,1);
            tmpl = sig(reftr, wb); %template: the most typical trial
            y = tmpl - mean(tmpl);
            for tr = 1:nresp
                if tr == reftr %the template trial stays where it is
                    trpk(tr) = 1;
                    continue
                end
                x = sig(tr,wb) - mean(sig(tr,wb)); %this trial, unshifted
                [r,lagv] = xcorr(x, y, maxlag, 'coeff'); %r peaks at +d when the trial is d late
                if all(isnan(r)) %a flat trial or template (no spikes), can't be aligned
                    continue
                end
                [pk,ip] = max(r);
                trlag(tr) = lagv(ip);
                trpk(tr) = pk;
            end

            %store, in ms
            if w == 1
                lags{1}{con} = trlag * bin_ms;
                lagpk{1}{con} = trpk;
                lagsd(1,con) = std(trlag * bin_ms);
            else
                lags_null{1}{con} = trlag * bin_ms;
                lagsd_null(1,con) = std(trlag * bin_ms);
            end
        end
    end
    clear sig Xw
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: Trial alignment results

%picture 1: how far each trial got shifted (typical-trial template), per condition, next to the same
%thing on the pre-stim null. Real structure = lags bunched up (small SD), noise = lags spread toward +/- max_lag_ms.
%picture 2 (one figure, a column per condition): the responsive-channel summed rate of every response trial,
%unaligned (top) and aligned (bottom), with the mean on top. A good alignment makes the mean sharper / taller.
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

    %picture 1: lag histograms, one per condition
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
    title(tl,[animal ' ' ProbeInfo.Areas{prb} '-Probe: trial lags, most typical trial as template (window ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms, max +/- ' num2str(max_lag_ms) ' ms)'],'Interpreter','none')

    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'align_lags']);
        % savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now
        exportgraphics(gcf,[save_directory '\Figures_StochRes\' animal '\' fname '.png'],'Resolution',png_res)
    end

    %picture 2: the trials unaligned (top row) and aligned (bottom row), one column per condition, all in one figure
    %(was one figure per condition)
    pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
    t_plot = t_bin(pidx); %ms
    pidx_ok = pidx(pidx > maxlag & pidx <= numel(t_bin) - maxlag); %bins that stay inside the data when shifted
    t_ok = t_bin(pidx_ok);
    panelnames = [{'unaligned'} align_methods];

    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(nmeth+1,ncond,'TileSpacing','compact','Padding','compact');
    for con = 1:ncond
        rtrs = find(resp{con} == true); %response trials
        nresp = numel(rtrs);
        if nresp < 2 %nothing to align: blank column with a note
            for p = 1:nmeth+1
                ax = nexttile(tl,(p-1)*ncond + con);
                axis(ax,'off')
                text(ax,0.5,0.5,[condis{con} ': fewer than 2 response trials'],'Units','normalized','HorizontalAlignment','center','Interpreter','none')
            end
            continue
        end
        sig = double(sum(rate{con}(rtrs,:,respchans),3)); % [response trials x all bins] spikes/s summed over the responsive channels
        ytop = max(sig(:,pidx),[],'all') * 1.05; %same y range for this condition's before and after panels

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

            ax = nexttile(tl,(p-1)*ncond + con); %row = version, column = condition
            plot(ax, t_ok, sh', 'Color',[.7 .7 .7], 'LineWidth',0.5) %every trial
            hold(ax,'on')
            plot(ax, t_ok, shmean, 'k', 'LineWidth',1.5) %their mean
            xline(ax,0,'r')
            xline(ax,align_use_ms(1),'b--')
            hold(ax,'off')
            xlim(ax,plot_ms)
            ylim(ax,[0 ytop])
            if p == nmeth+1
                xlabel(ax,'ms from stim onset')
            end
            if con == 1
                ylabel(ax,{panelnames{p}, 'spikes/s, summed resp. ch.'},'Interpreter','none')
            end
            %how wiggly the mean is in the alignment window (SD over time): goes up when trials line up better
            title(ax,{[condis{con} ', ' panelnames{p}], ['mean''s SD in window: ' num2str(std(shmean(aidx)),'%.1f') ' spikes/s']},'Interpreter','none')
        end
    end
    title(tl,[animal ' ' ProbeInfo.Areas{prb} '-Probe: response trials (grey) and mean (black), before (top) and after (bottom) alignment (blue dashed = start of the alignment window)'],'Interpreter','none')

    if plotsave == 1
        % fname = sprintf([animal '-' condis{con} '-' ProbeInfo.Areas{prb} '_Probe' '-' 'align_trials']); %commented out: one figure per condition before
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'align_trials']);
        % savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now
        exportgraphics(gcf,[save_directory '\Figures_StochRes\' animal '\' fname '.png'],'Resolution',png_res)
    end
    clear sig sh
end




%Now beginning the aligned-trial analyses

%% Data Analysis: Trial-to-trial correlation on aligned trials

%Redo the trial consistency with the trials lined up (typical-trial template from the alignment section),
%next to the same thing unaligned, in the alignment window (align_use_ms).
%(aligned evoked / induced spectra were tried and dropped)
%Trial-to-trial correlation (as NNx_LFPCorrelation_LvsLW.m): every pair of response trials, their smoothed
%rate time courses correlated (r), Fisher z, averaged within each condition, back to r (a group's
%conditions are never pooled). Group 2 vs group 1, per condition pair (e.g. L_8 vs LW_8): the labels of
%that pair's trials are shuffled n_perm times, p = how often a shuffle gives a difference at least as big
%as the real one.
%Null: the same consistency numbers on the same-length pre-stim stretch, aligned with the null lags, so
%you can see how much "consistency" turns up with no stimulus (lining trials up always adds some).
if do_align == 1
    disp('MUA: Trial consistency, unaligned vs aligned')
    tic

    al_names = {'unaligned','typical'}; %the versions compared
    al_meth = [0 find(strcmp(align_methods,'typical'))]; %which alignment method's lags each version uses (0 = no shift)
    nver = numel(al_names);
    wbr = winbins{1}; %bins of the alignment window (from the alignment section)
    wbn = winbins{2}; %bins of its pre-stim copy (null)
    nwin_al = numel(wbr);

    tc_trial = cell(1,nver); %{versions}{conditions} [response trials x channels] each trial's mean r with the other trials of its condition
    % tc_grp = cell(1,nver); %commented out: pooled each group's conditions, conditions stay apart now
    % tc_grp_null = cell(1,nver);
    % ngrp = numel(condis_grp);
    tc_cond = cell(1,nver); %{versions} [channels x conditions] mean pairwise trial r within each condition
    tc_cond_null = cell(1,nver); %{versions} [channels x conditions] the same on the pre-stim window
    tc_p = cell(1,nver); %{versions} [channels x pairs] permutation p, group 2 vs group 1 condition of each pair

    for v = 1:nver
        disp(['MUA: ' al_names{v}])
        tc_trial{v} = cell(1,ncond);

        %line every condition's response trials up (or not) and stack them
        allr = []; % [all response trials x window bins x channels] smoothed rate, real window
        allrn = []; % [all response trials x window bins x channels] smoothed rate, pre-stim window
        condlab = []; %condition of each stacked trial
        for con = 1:ncond
            rtrs = find(resp{con} == true); %response trials
            nresp = numel(rtrs);
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
            r_aln = zeros(nresp,nwin_al,nch,'single'); %rate, pre-stim window
            for tr = 1:nresp
                r_al(tr,:,:) = rate{con}(rtrs(tr), wbr + lb(tr), :);
                r_aln(tr,:,:) = rate{con}(rtrs(tr), wbn + lbn(tr), :);
            end
            allr = cat(1, allr, r_al);
            allrn = cat(1, allrn, r_aln);
            condlab = [condlab; con*ones(nresp,1)];
        end
        % grplab = condgrp(condlab)'; %commented out: groups aren't pooled any more
        ntot = numel(condlab);

        %trial-to-trial correlation, per channel
        tc_cond{v} = nan(nch,ncond);
        tc_cond_null{v} = nan(nch,ncond);
        tc_p{v} = nan(nch,npair);
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

                %each trial with the other trials of its condition (real window)
                if w == 1
                    for con = 1:ncond
                        idx = condlab == con;
                        if sum(idx) < 2
                            continue
                        end
                        tc_trial{v}{con}(:,ch) = tanh(mean(Z(idx,idx),2,'omitnan'));
                    end
                end

                %per condition (a group's conditions are never pooled)
                zc = nan(1,ncond);
                for con = 1:ncond
                    ci = condlab == con;
                    if sum(ci) >= 2
                        zc(con) = mean(Z(ci,ci),'all','omitnan');
                    end
                end
                if w == 1
                    tc_cond{v}(ch,:) = tanh(zc);
                else
                    tc_cond_null{v}(ch,:) = tanh(zc);
                end

                %group 2 vs group 1 inside each condition pair, shuffling only that pair's trial labels (real window only)
                if w == 1
                    for k = 1:npair
                        ca = pairs(k,1); %group 1 condition
                        cb = pairs(k,2); %group 2 condition
                        if isnan(zc(ca)) || isnan(zc(cb))
                            continue
                        end
                        pti = find(condlab == ca | condlab == cb); %the pair's trials in the stack
                        Zp = Z(pti,pti);
                        plab = condlab(pti); %their condition
                        obs = zc(cb) - zc(ca);
                        % nulldiff = zeros(n_perm,1); %commented out (with the loop below): all shuffles in one matrix product now, ~5x faster, same numbers
                        % for p = 1:n_perm
                        %     shlab = plab(randperm(numel(plab))); %shuffled labels within the pair
                        %     ia = shlab == ca;
                        %     ib = shlab == cb;
                        %     nulldiff(p) = mean(Zp(ib,ib),'all','omitnan') - mean(Zp(ia,ia),'all','omitnan');
                        % end
                        %the mean of Zp over a set of trials = (sum over the set's pairs) / (number of those pairs that
                        %aren't NaN), and for every shuffle at once that is one matrix product (same random draws as before)
                        ntp = numel(plab);
                        IB = zeros(ntp,n_perm); % [pair trials x shuffles] 1 = labelled as condition B in that shuffle
                        for p = 1:n_perm
                            shlab = plab(randperm(ntp)); %shuffled labels within the pair
                            IB(:,p) = shlab == cb;
                        end
                        IA = 1 - IB; %labelled as condition A
                        Wp = double(isnan(Zp) == false); %1 = a pair that counts
                        Z0 = Zp;
                        Z0(Wp == 0) = 0;
                        nulldiff = (sum(IB .* (Z0*IB),1) ./ sum(IB .* (Wp*IB),1) - sum(IA .* (Z0*IA),1) ./ sum(IA .* (Wp*IA),1))'; % [shuffles x 1]
                        tc_p{v}(ch,k) = (sum(abs(nulldiff) >= abs(obs)) + 1) / (n_perm + 1);
                    end
                end
            end
        end
    end
    clear allr allrn r_al r_aln X C Z Zp Z0 Wp IA IB
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: Aligned-trial results

%per condition pair (slider; a group's conditions are never pooled):
%picture 1: trial consistency on the probe (one square per channel): condition A (group 1) r, condition B
%           (group 2) r, B - A (black square = p < 0.05), and the same difference on the pre-stim null; unaligned and aligned
%picture 2: how much lining the trials up raised the consistency, per condition, next to the same gain on the null
%picture 3: each trial's consistency with its condition (mean over the responsive channels), unaligned vs aligned,
%           conditions side by side in their pairs
if do_align == 1 && plot_corr == 1
    disp('MUA: Making aligned-trial pictures')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    [nrow,nshank] = size(pmap);
    vun = find(strcmp(al_names,'unaligned'));
    vstr = find(strcmp(al_names,'typical'));
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
        %pictures 1 and 2: probe maps, ONE figure each with a slider over the condition pairs (a panel per
        %pair, its group 1 condition A vs its group 2 condition B; the groups' conditions are never pooled).
        %Each panel is a list of maps: values (channels x 1), title, diverging colours or not, channels to
        %mark with a black square. Colour ranges are shared by all pairs, so the panels compare
        for fg = 1:2
            allpv = cell(1,npair); %{pairs}{maps} values
            allpt = cell(1,npair); %{pairs}{maps} titles
            allpm = cell(1,npair); %{pairs}{maps} channels to mark
            for k = 1:npair
                ca = pairs(k,1); %group 1 condition of this pair
                cb = pairs(k,2); %group 2 condition
                pv = {};
                pt = {};
                pd = [];
                pm = {};
                if fg == 1
                    for v = 1:nver
                        pv(end+1:end+4) = {tc_cond{v}(:,ca), tc_cond{v}(:,cb), tc_cond{v}(:,cb) - tc_cond{v}(:,ca), tc_cond_null{v}(:,cb) - tc_cond_null{v}(:,ca)};
                        pt(end+1:end+4) = {[al_names{v} ': ' condis{ca} ' trial r'], [al_names{v} ': ' condis{cb} ' trial r'], ...
                            [al_names{v} ': ' condis{cb} ' - ' condis{ca} ' (square = p < 0.05)'], [al_names{v} ': ' condis{cb} ' - ' condis{ca} ', pre-stim null']};
                        pd(end+1:end+4) = [0 0 1 1];
                        pm(end+1:end+4) = {[], [], tc_p{v}(:,k) < 0.05, []};
                    end
                else
                    pv = {tc_cond{vstr}(:,ca) - tc_cond{vun}(:,ca), tc_cond{vstr}(:,cb) - tc_cond{vun}(:,cb), ...
                        tc_cond_null{vstr}(:,ca) - tc_cond_null{vun}(:,ca), tc_cond_null{vstr}(:,cb) - tc_cond_null{vun}(:,cb)};
                    pt = {[condis{ca} ': gain from aligning (aligned - unaligned r)'], [condis{cb} ': gain from aligning (aligned - unaligned r)'], ...
                        [condis{ca} ': same gain on the pre-stim null'], [condis{cb} ': same gain on the pre-stim null']};
                    pd = [1 1 1 1];
                    pm = {[], [], [], []};
                end
                allpv{k} = pv;
                allpt{k} = pt;
                allpm{k} = pm;
            end
            if fg == 1
                tlsize = [nver 4];
                figtitle = [animal ' ' ProbeInfo.Areas{prb} '-Probe: trial-to-trial consistency (mean pairwise r of the smoothed rate, ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms)'];
                ftag = 'tc_maps';
            else
                tlsize = [1 4];
                figtitle = [animal ' ' ProbeInfo.Areas{prb} '-Probe: how much lining the trials up raised their consistency (real vs no-stim null)'];
                ftag = 'tc_gain';
            end
            % (the old per-group version, tc_grp group 1 vs group 2 in one figure, was dropped: it pooled conditions)

            %colour ranges over every pair: one for the r panels, one symmetric one for the difference panels
            rvals = [];
            dvals = [];
            for k = 1:npair
                rvals = [rvals; cell2mat(allpv{k}(pd == 0)')];
                dvals = [dvals; cell2mat(allpv{k}(pd == 1)')];
            end
            rlim = [min(rvals,[],'omitnan') max(rvals,[],'omitnan')];
            if isempty(rlim) || any(isnan(rlim)) || rlim(1) == rlim(2)
                rlim = [0 1];
            end
            dmax = max(abs(dvals),[],'omitnan');
            if isempty(dmax) || isnan(dmax) || dmax == 0
                dmax = 0.01;
            end

            mapfig = figure();
            set(mapfig, 'Position', get(0, 'Screensize'));
            for k = 1:npair
                %this pair's panel (only the first starts visible); Tag + UserData let the slider find it
                pan = uipanel(mapfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
                    'Tag','condpanel', 'UserData',k, 'Visible', k == 1);
                tl = tiledlayout(pan,tlsize(1),tlsize(2),'TileSpacing','compact','Padding','compact');
                pv = allpv{k};
                pt = allpt{k};
                pm = allpm{k};
                for kk = 1:numel(pv)
                    ax = nexttile(tl,kk);
                    zmap = nan(size(pmap)); % [rows x shanks]
                    for ch = 1:nch
                        zmap(chrow(ch),chshank(ch)) = pv{kk}(ch);
                    end
                    h = imagesc(ax, zmap);
                    set(h,'AlphaData',isnan(zmap) == false) %no site / removed = blank
                    if pd(kk) == 1
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
                    if isempty(pm{kk}) == false && any(pm{kk})
                        plot(ax, chshank(pm{kk}), chrow(pm{kk}), 'ks', 'MarkerSize',16, 'LineWidth',1.5) %p < 0.05
                    end
                    for i = 1:numel(rrow)
                        text(ax, rshank(i), rrow(i), 'x', 'HorizontalAlignment','center','Color',[.5 .5 .5])
                    end
                    hold(ax,'off')
                    set(ax,'XTick',1:nshank,'YTick',1:nrow)
                    xlabel(ax,'shank')
                    ylabel(ax,'row (1 = top)')
                    title(ax, pt{kk}, 'Interpreter','none')
                end
                title(tl,{figtitle, [pairnames{k} ' (pair ' num2str(k) ' of ' num2str(npair) '); bold = responsive channel, x = removed']},'Interpreter','none')
            end

            %the slider (finds the panels by their Tag, so it still works in the saved figure)
            if npair > 1
                uicontrol(mapfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
                    'String','pair:', 'HorizontalAlignment','right')
                uicontrol(mapfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
                    'Min',1, 'Max',npair, 'Value',1, 'SliderStep',[1 1]/(npair-1), ...
                    'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
                    num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
            end
            drawnow %show it now; the script doesn't wait, it carries on

            if plotsave == 1
                fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' ftag]);
                % savefig(mapfig,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now, one per panel
                %one PNG per panel (= pair): show only panel k, save its layout (not the slider), then back to panel 1
                ptags = strrep(pairnames,' ','_'); %file name tag of each panel, e.g. L_8_vs_LW_8
                pans = findobj(mapfig,'Tag','condpanel');
                for k = 1:numel(ptags)
                    set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
                    pk = findobj(pans,'flat','UserData',k);
                    tlk = findobj(pk,'Type','tiledlayout');
                    exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' ptags{k} '.png'],'Resolution',png_res)
                end
                set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
                set(findobj(mapfig,'Style','slider'),'Value',1);
            end
        end

        %picture 3: each trial's consistency with its condition, averaged over the responsive channels
        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        ax = axes();
        hold(ax,'on')
        xt = [];
        xtl = {};
        porder = reshape(pairs',1,[]); %conditions in pair order: pair 1 A, pair 1 B, pair 2 A, ...
        for j = 1:ncond
            con = porder(j);
            trval = cell(1,nver); %each trial's value per version
            for v = 1:nver
                trval{v} = mean(tc_trial{v}{con}(:,respchans),2,'omitnan'); % [response trials x 1]
            end
            % xpos = (con-1)*3 + [1 2]; %commented out: conditions now sit in pairs, with a gap between pairs
            xpos = (j-1)*3 + floor((j-1)/2)*2 + [1 2]; %unaligned, aligned
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
        title(ax,{[animal ' ' ProbeInfo.Areas{prb} '-Probe: per-trial consistency, unaligned vs aligned (typical-trial template)'], ...
            ['conditions side by side in their pairs (' strjoin(pairnames,', ') '); grey lines join the same trial before and after aligning']},'Interpreter','none')

        if plotsave == 1
            fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'tc_trials']);
            % savefig(gcf,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now
            exportgraphics(gcf,[save_directory '\Figures_StochRes\' animal '\' fname '.png'],'Resolution',png_res)
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

    g2cons = find(condgrp == 2); %group 2's conditions, each analysed on its own (never pooled)
    rch = find(respchans); %responsive V1 channels
    nrch = numel(rch);
    nch_s1 = numel(chanlist_s1);
    if exist('bwconncomp','file') == 0
        error('The shape-map clusters need bwconncomp (Image Processing Toolbox)')
    end

    %sliding blocks for the shape match (see below) and the cluster threshold, the same for every condition
    hb = round(shape_block_ms/2/tf_step_ms); %windows either side of the block centre
    L = 2*hb + 1; %windows per block
    bc = 1+hb : nt-hb; %block centres (window index)
    nbk = numel(bc);
    zthr = norminv(1 - shuf_p/2); %pixel |z| threshold for the clusters

    %results per condition {1 x conditions}, empty / NaN for group 1 (only group 2 is analysed)
    % ntr2 = 0; %commented out: group 2's conditions used to be pooled into one analysis
    ntr_s1 = zeros(1,ncond); %response trials used per condition
    v1tr = cell(1,ncond); %{conditions} [freqs x windows x trials] V1 rhythmicity, dB
    r_shape_ch = cell(1,ncond); %{conditions} [freqs x windows x S1 channels] shape r
    shape_zpix = cell(1,ncond); %{conditions} [freqs x windows x S1 channels] pixel z vs the shuffles
    shape_ppix = cell(1,ncond); %{conditions} [freqs x windows x S1 channels] pixel p
    shape_sig = cell(1,ncond); %{conditions} [freqs x windows x S1 channels] pixels in significant clusters
    shape_clustthr = nan(1,ncond); %cluster threshold per condition

    for con = g2cons
        nresp = size(lfp_s1{con},1); %response trials (lfp_s1 only holds these)
        if nresp < 3 %the trial shuffles need at least 3 trials
            run_notes(end+1) = string(condis{con}) + ": fewer than 3 response trials, no S1 shape match";
            continue
        end
        rtrs = find(resp{con} == true); %their rows in spkcount
        ntr_s1(con) = nresp;
        disp(['S1: ' condis{con} ' (' num2str(nresp) ' response trials)'])

        %per-trial values, this condition only
        v1trsum = zeros(nftf,nt,nresp); % [freqs x windows x trials] V1 rhythmicity summed over the responsive channels
        v1trn = zeros(nftf,nt,nresp); %how many responsive channels went into each
        s1tr = zeros(nftf,nt,nresp,nch_s1,'single'); % [freqs x windows x trials x S1 channels] S1 power, single saves memory

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
                v1trsum(:,ti,:) = v1trsum(:,ti,:) + reshape(R, nftf, 1, nresp); %per trial
                v1trn(:,ti,:) = v1trn(:,ti,:) + reshape(good, nftf, 1, nresp);
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
                s1tr(:,ti,:,c1) = single(reshape(P(fk_tf,:), nftf, 1, nresp)); %per trial
            end
        end

        %V1 rhythmicity per trial: responsive-channel mean, dB above random firing (NaN = no spikes)
        v1c = 10*log10(v1trsum ./ v1trn); % [freqs x windows x trials]
        v1c(isinf(v1c)) = NaN;
        v1tr{con} = v1c;

        %shape of the time courses: at each frequency, within each trial, correlate the V1 rhythmicity curve and
        %the S1 power curve over a sliding block of windows (shape_block_ms of centres), then average the trials'
        %r (Fisher z). Pearson ignores offset and scale, so a curve sitting below its pre-stim level (negative dB)
        %still counts if it rises and falls with the other one. Neighbouring windows overlap a lot, so each block
        %holds only a few independent points: r is noisy and large values happen by chance -> trial-shuffle null.
        %
        %Trial-shuffle null: V1 curves of each trial paired with ANOTHER trial's S1 curves (no trial keeps its own
        %partner; trials of this condition only), whole map recomputed, n_shuf times. Whatever survives shuffling
        %(the time course all trials share, the overlapping windows) is chance level. Each pixel gets a z against
        %its shuffles; neighbouring pixels with |z| past the shuf_p threshold form clusters (positive and negative
        %apart), cluster mass = sum of |z|; a real cluster counts when its mass beats the (1 - clust_p) percentile
        %of the biggest chance cluster of each shuffle, over ALL S1 channels (so the whole probe is covered).
        %Each condition gets its own null and threshold.
        %Speed: each trial's block curves are centred and scaled once (for V1 and for S1), so a pairing is just
        %a multiply and a sum. Windows where V1 had no spikes are left out of the V1 curve but not the S1 curve
        %(slight approximation; rare once the responsive channels are averaged)

        %V1: each trial's block curve, centred and scaled to length 1 (once, the same for every S1 channel)
        Vn = zeros(nftf,L,nresp,nbk,'single'); % [freqs x block windows x trials x blocks]
        Vok = false(nftf,1,nresp,nbk); %trial usable in that block (>= 5 windows with spikes, not flat)
        for k = 1:nbk
            a = v1c(:, bc(k)-hb : bc(k)+hb, :); % [freqs x block windows x trials]
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
        rng('default') %same shuffles every run (per condition, so each is reproducible on its own)
        permset = zeros(n_shuf, nresp);
        for s = 1:n_shuf
            pm = randperm(nresp);
            while any(pm == 1:nresp) %no trial paired with itself
                pm = randperm(nresp);
            end
            permset(s,:) = pm;
        end

        %Speed: every V1 trial x S1 trial pairing is worked out ONCE per S1 channel (one batched matrix product,
        %pagemtimes), so the real map is the trials matched with themselves and each shuffle just picks other
        %entries of that table, instead of redoing the multiply-and-sum for every shuffle (~4x faster, same numbers)
        V2 = reshape(permute(Vn,[3 2 1 4]), nresp, L, nftf*nbk); % [trials x block windows x (freqs*blocks)] V1 curves
        Vokp = reshape(permute(Vok,[3 1 4 2]), nresp, nftf*nbk); % [trials x (freqs*blocks)] usable
        diagidx = (1:nresp)' + ((1:nresp)'-1)*nresp; %table entries of each trial with itself
        shufidx = (1:nresp) + (permset-1)*nresp; % [shuffles x trials] table entries of trial i with its shuffled partner

        %worker pool (Parallel Computing Toolbox only; otherwise nw = 0 = one core, a normal loop in order)
        nw = 0;
        if n_workers > 0 && license('test','Distrib_Computing_Toolbox') && isempty(ver('parallel')) == false
            pool = gcp('nocreate'); %use the pool if one is open
            if isempty(pool)
                pool = parpool; %default size (first time in a session takes a while)
            end
            nw = min(n_workers, pool.NumWorkers);
        end
        cname = condis{con}; %this condition's name, for the progress lines

        rsh = nan(nftf,nt,nch_s1); % [freqs x windows (block centres) x S1 channels], NaN at the edges
        zpx = nan(nftf,nt,nch_s1); % [freqs x windows x S1 channels] pixel z vs the shuffles
        ppx = nan(nftf,nt,nch_s1); % [freqs x windows x S1 channels] pixel p (two-sided, from the shuffles)
        nullmax = zeros(n_shuf, nch_s1); %biggest chance cluster mass per shuffle and channel
        % for c1 = 1:nch_s1 %commented out: now on nw workers (one S1 channel each; nw = 0 runs it as this loop)
        parfor (c1 = 1:nch_s1, nw) %S1 channels don't depend on each other (each worker gets V2 and the index tables, a few MB)
            s1l = 10*log10(double(s1tr(:,:,:,c1))); % [freqs x windows x trials] dB
            s1l(isinf(s1l)) = NaN;

            %S1: each trial's block curve, centred and scaled to length 1
            Sn = zeros(nftf,L,nresp,nbk,'single');
            for k = 1:nbk
                b = s1l(:, bc(k)-hb : bc(k)+hb, :);
                b(isnan(b)) = 0;
                db = b - mean(b,2);
                db = db ./ sqrt(sum(db.^2,2));
                db(isnan(db)) = 0;
                Sn(:,:,:,k) = single(db);
            end

            %the table: r of V1 trial i with S1 trial j, every pairing, for every frequency and block, Fisher z
            S2 = reshape(permute(Sn,[2 3 1 4]), L, nresp, nftf*nbk); % [block windows x trials x (freqs*blocks)]
            A = double(reshape(pagemtimes(V2, S2), nresp*nresp, nftf*nbk)); % [(trials*trials) x (freqs*blocks)]
            A(A > 0.9999) = 0.9999; %keep atanh finite
            A(A < -0.9999) = -0.9999;
            A = atanh(A);

            %real map: each trial's V1 with its own S1, Fisher-z mean over trials
            % rt = double(sum(Vn .* Sn, 2)); %commented out (and the 3 lines after): the table above does this in one go
            % rt(Vok == false) = NaN;
            % rt(rt > 0.9999) = 0.9999;
            % rt(rt < -0.9999) = -0.9999;
            % zreal = reshape(mean(atanh(rt),3,'omitnan'), nftf, nbk);
            rt = A(diagidx,:); % [trials x (freqs*blocks)]
            rt(Vokp == false) = NaN;
            zreal = reshape(mean(rt,1,'omitnan'), nftf, nbk); % [freqs x blocks]

            %shuffled maps
            zsh = zeros(nftf,nbk,n_shuf,'single');
            for s = 1:n_shuf
                % rt = double(sum(Vn .* Sn(:,:,permset(s,:),:), 2)); %commented out: picked from the table instead
                rt = A(shufidx(s,:),:); %V1 of trial i with S1 of another trial
                rt(Vokp == false) = NaN;
                zsh(:,:,s) = single(reshape(mean(rt,1,'omitnan'), nftf, nbk));
            end
            zsh(linemask_tf,:,:) = NaN;
            zreal(linemask_tf,:) = NaN;

            %pixel z and p against the shuffles
            mu = mean(zsh,3,'omitnan');
            sd = std(zsh,0,3,'omitnan');
            zp = (zreal - mu) ./ sd; % [freqs x blocks]
            pp = (sum(abs(zsh - mu) >= abs(zreal - mu), 3) + 1) / (n_shuf + 1);

            %biggest chance cluster of each shuffle (its own map against the same shuffle mean / SD)
            nm = zeros(n_shuf,1);
            for s = 1:n_shuf
                zs = (double(zsh(:,:,s)) - mu) ./ sd;
                mx = 0;
                for sgn = [1 -1] %positive and negative clusters apart
                    cc = bwconncomp(sgn*zs > zthr, 4);
                    for q = 1:cc.NumObjects
                        mx = max(mx, sum(abs(zs(cc.PixelIdxList{q}))));
                    end
                end
                nm(s) = mx;
            end

            %this channel's results first, then into the big arrays in one go (parfor needs it that way)
            rsh_c = nan(nftf,nt);
            zpx_c = nan(nftf,nt);
            ppx_c = nan(nftf,nt);
            rsh_c(:,bc) = tanh(zreal);
            zpx_c(:,bc) = zp;
            ppx_c(:,bc) = pp;
            rsh(:,:,c1) = rsh_c;
            zpx(:,:,c1) = zpx_c;
            ppx(:,:,c1) = ppx_c;
            nullmax(:,c1) = nm;
            disp(['S1: ' cname ' shape-map shuffles done for CH' num2str(chanlist_s1(c1))])
        end

        %cluster threshold: biggest chance cluster over all S1 channels, per shuffle
        cthr = prctile(max(nullmax,[],2), 100*(1 - clust_p));
        sgx = false(nftf,nt,nch_s1); %pixels in significant clusters
        for c1 = 1:nch_s1
            zp = zpx(:,bc,c1);
            sigm = false(nftf,nbk);
            for sgn = [1 -1]
                cc = bwconncomp(sgn*zp > zthr, 4);
                for q = 1:cc.NumObjects
                    if sum(abs(zp(cc.PixelIdxList{q}))) > cthr
                        sigm(cc.PixelIdxList{q}) = true;
                    end
                end
            end
            sgx(:,bc,c1) = sigm;
        end

        %store this condition's results
        r_shape_ch{con} = rsh;
        shape_zpix{con} = zpx;
        shape_ppix{con} = ppx;
        shape_sig{con} = sgx;
        shape_clustthr(con) = cthr;
    end
    if all(ntr_s1 == 0)
        error(['No ' grpnames{2} ' condition has 3 or more response trials, nothing to compare'])
    end
    clear X P Y R v1trsum v1trn s1tr s1l a b da db rt Vn Sn zsh okv v1c rsh zpx ppx sgx V2 Vokp shufidx
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
    %over a shape_block_ms block, Fisher-z averaged over the condition's trials).
    %ONE figure, a panel per group 2 condition (slider; conditions are never pooled), one shared colour range
    %(an across-trial r per pixel was tried first and dropped: it compared values, not the shape over time)
    s1cons = find(ntr_s1 > 0); %group 2 conditions that were analysed
    rlim = 0;
    for con = s1cons
        rlim = max(rlim, max(abs(r_shape_ch{con}(:)),[],'omitnan')); %one symmetric colour range for every panel
    end
    if isnan(rlim) || rlim == 0
        rlim = 1;
    end

    s1fig = figure();
    set(s1fig, 'Position', get(0, 'Screensize'));
    for k = 1:numel(s1cons)
        con = s1cons(k);
        rmap = r_shape_ch{con};
        sigmap = shape_sig{con};
        %this condition's panel (only the first starts visible); Tag + UserData let the slider find it
        pan = uipanel(s1fig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
            'Tag','condpanel', 'UserData',k, 'Visible', k == 1);
        tl = tiledlayout(pan,nrow_s1,nshank_s1,'TileSpacing','none','Padding','tight');
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
                set(h,'AlphaData',(isnan(rmap(:,:,c1)) == false) .* (0.35 + 0.65*sigmap(:,:,c1))) %NaN = blank, not significant = faded
                axis(ax,'xy')
                colormap(ax,divmap)
                clim(ax,[-rlim rlim])
                hold(ax,'on')
                xline(ax,0,'k')
                if any(sigmap(:,:,c1),'all')
                    contour(ax, tf_t, tf_f, double(sigmap(:,:,c1)), [0.5 0.5], 'k', 'LineWidth',1) %outline the significant clusters
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
        % title(tl,{[animal ' S1 probe, ' grpnames{2} ' (' num2str(ntr2) ' response trials): ...  %commented out: one pooled group 2 map before
        title(tl,{[animal ' S1 probe, ' condis{con} ' (' num2str(ntr_s1(con)) ' response trials; condition ' num2str(k) ' of ' num2str(numel(s1cons)) '): shape match between V1 spike rhythmicity (responsive-channel mean) and each S1 channel''s power over time, per frequency'], ...
            ['r over a ' num2str(shape_block_ms) ' ms block of ' num2str(tf_win_ms) ' ms windows around each time, per trial, Fisher-z averaged; red = they rise and fall together; bold blue = CH' num2str(s1_refch)], ...
            ['outlined / not faded = cluster beats ' num2str(n_shuf) ' trial shuffles of this condition (pixel p < ' num2str(shuf_p) ', cluster p < ' num2str(clust_p) ', corrected over the whole S1 probe)']},'Interpreter','none')
    end

    %the slider (finds the panels by their Tag, so it still works in the saved figure)
    if numel(s1cons) > 1
        uicontrol(s1fig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
            'String','condition:', 'HorizontalAlignment','right')
        uicontrol(s1fig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
            'Min',1, 'Max',numel(s1cons), 'Value',1, 'SliderStep',[1 1]/(numel(s1cons)-1), ...
            'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
            num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
    end
    drawnow %show it now; the script doesn't wait, it carries on

    if plotsave == 1
        fname = sprintf([animal '-' 'S1_vs_V1_timefreq_rmaps_S1probe']);
        % savefig(s1fig,[save_directory '\Figures_StochRes\' animal '\' fname],'compact') %commented out: PNGs for now, one per panel
        %one PNG per panel (= group 2 condition): show only panel k, save its layout (not the slider), then back to panel 1
        ptags = condis(s1cons); %file name tag of each panel
        pans = findobj(s1fig,'Tag','condpanel');
        for k = 1:numel(ptags)
            set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
            pk = findobj(pans,'flat','UserData',k);
            tlk = findobj(pk,'Type','tiledlayout');
            exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' ptags{k} '.png'],'Resolution',png_res)
        end
        set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
        set(findobj(s1fig,'Style','slider'),'Value',1);
    end
end




%Now beginning the wrap-up

%% Post-run report

%everything worth knowing about this run in one place, at the end (nothing was printed along the way).
%It's saved with the results too, as `report`. Parts that were switched off are just skipped
report = strings(0,1);
report(end+1) = "===== StochRes MUA report: " + animal + " (" + string(datetime('now','Format','yyyy-MM-dd HH:mm')) + ") =====";
report(end+1) = "Probe: " + ProbeInfo.Areas{prb} + " (probe " + prb + "), " + nch + " kept channels";
for g = 1:numel(condis_grp)
    report(end+1) = "Group " + g + " (" + grpnames{g} + "): " + strjoin(string(condis_grp{g}), ", ");
end
report(end+1) = "Files:";
for i = 1:size(filesused,1)
    report(end+1) = "  " + filesused(i,1) + ": " + filesused(i,2) + " (" + filesused(i,3) + " good trials, " + filesused(i,4) + " left out by tr_conditional_use)";
end

%response trials
if exist('resp','var')
    report(end+1) = "Response trials (template matching, " + resp_ms(1) + "-" + resp_ms(2) + " ms; score > the " + resp_null_prct + "th percentile of the condition's fake no-stim scores):";
    for con = 1:ncond
        report(end+1) = "  " + condis{con} + ": " + sum(resp{con}) + " of " + ntrials(con) + " (" + sprintf('%.0f',100*respfrac(con)) + "%), threshold " + sprintf('%.2f',respthr(con)) + ...
            ", " + sum(respchans_cond(:,con)) + " responsive channels, " + tmpl_rounds(con) + " template round(s)";
    end
    % for g = 1:numel(condis_grp) %commented out: pooled each group's conditions; now compared per pair
    %     report(end+1) = "  group " + grpnames{g} + ": " + sprintf('%.0f',100*respfrac_grp(g)) + "%";
    % end
    for k = 1:npair
        ca = pairs(k,1);
        cb = pairs(k,2);
        report(end+1) = "  " + pairnames{k} + ": " + sum(resp{ca}) + " of " + ntrials(ca) + " (" + sprintf('%.0f',100*respfrac(ca)) + "%) vs " + ...
            sum(resp{cb}) + " of " + ntrials(cb) + " (" + sprintf('%.0f',100*respfrac(cb)) + "%)";
    end
    report(end+1) = "  responsive channels, any condition (" + sum(respchans) + "): " + strjoin(string(chanlist(respchans)'), ", ");

    %Quality Control details (keep_qc): every condition's trials by call, and its own responsive channels
    if keep_qc == 1
        report(end+1) = "Quality Control (trials as file number within the condition : trial number in that file):";
        for con = 1:ncond
            trlab = "f" + string(trinfo{con}(:,1)) + ":" + string(trinfo{con}(:,2)); % [trials x 1] label of each good trial
            report(end+1) = "  " + condis{con} + " kept (response, " + sum(resp{con}) + "): " + strjoin(trlab(resp{con})', ", ");
            report(end+1) = "  " + condis{con} + " left out (no response, " + sum(resp{con} == false) + "): " + strjoin(trlab(resp{con} == false)', ", ");
            report(end+1) = "  " + condis{con} + " responsive channels (" + sum(respchans_cond(:,con)) + "): " + strjoin(string(chanlist(respchans_cond(:,con))'), ", ");
        end
    end
end

%alignment jitter
if do_align == 1 && exist('lagsd','var')
    report(end+1) = "Alignment lag jitter, SD in ms (pre-stim null in brackets), window " + align_use_ms(1) + "-" + align_use_ms(2) + " ms:";
    for m = 1:numel(align_methods)
        txt = "  " + align_methods{m} + ":";
        for con = 1:ncond
            txt = txt + "  " + condis{con} + " " + sprintf('%.1f',lagsd(m,con)) + " (" + sprintf('%.1f',lagsd_null(m,con)) + ")";
        end
        report(end+1) = txt;
    end
    %Quality Control (keep_qc): which trial each condition was aligned to
    if keep_qc == 1
        txt = "  reference (most typical) trial, file:trial (r with the others' mean):";
        for con = 1:ncond
            if isnan(align_ref(con))
                txt = txt + "  " + condis{con} + " none";
            else
                txt = txt + "  " + condis{con} + " f" + trinfo{con}(align_ref(con),1) + ":" + trinfo{con}(align_ref(con),2) + " (" + sprintf('%.2f',align_ref_r(con)) + ")";
            end
        end
        report(end+1) = txt;
    end
end

%trial-to-trial consistency
% if do_align == 1 && exist('tc_grp','var') %commented out: the per-group version pooled conditions
if do_align == 1 && exist('tc_cond','var')
    report(end+1) = "Trial-to-trial consistency per condition pair, mean pairwise r over the responsive channels (pre-stim null in brackets):";
    for k = 1:npair
        ca = pairs(k,1);
        cb = pairs(k,2);
        for v = 1:numel(al_names)
            txt = "  " + pairnames{k} + ", " + al_names{v} + ":  " + ...
                condis{ca} + " " + sprintf('%.3f',mean(tc_cond{v}(respchans,ca),'omitnan')) + " (" + sprintf('%.3f',mean(tc_cond_null{v}(respchans,ca),'omitnan')) + ")  " + ...
                condis{cb} + " " + sprintf('%.3f',mean(tc_cond{v}(respchans,cb),'omitnan')) + " (" + sprintf('%.3f',mean(tc_cond_null{v}(respchans,cb),'omitnan')) + ")" + ...
                ";  channels with p < 0.05: " + sum(tc_p{v}(:,k) < 0.05) + " of " + nch;
            report(end+1) = txt;
        end
    end
end

%S1 vs V1 shape match, per group 2 condition
if do_S1 == 1 && exist('shape_sig','var')
    report(end+1) = "S1 vs V1 shape match (" + n_shuf + " shuffles per condition; single S1 channel CH" + s1_refch + "):";
    for con = find(ntr_s1 > 0)
        report(end+1) = "  " + condis{con} + " (" + ntr_s1(con) + " response trials, cluster threshold " + sprintf('%.1f',shape_clustthr(con)) + "):";
        nsig = 0;
        for c1 = 1:numel(chanlist_s1)
            sigc = shape_sig{con}(:,:,c1);
            if any(sigc,'all')
                nsig = nsig + 1;
                rr = r_shape_ch{con}(:,:,c1);
                rr(sigc == false) = NaN; %only the significant pixels
                [~,ip] = max(abs(rr(:))); %strongest one (max skips NaN)
                [fi,ti] = ind2sub(size(rr), ip);
                report(end+1) = "    CH" + chanlist_s1(c1) + ": significant, peak r " + sprintf('%.2f',rr(ip)) + " at " + sprintf('%.0f',tf_f(fi)) + " Hz, " + tf_t(ti) + " ms";
            end
        end
        report(end+1) = "    " + nsig + " of " + numel(chanlist_s1) + " S1 channels with a significant cluster";
    end
end

%notes collected along the way
report(end+1) = "Notes:";
if isempty(run_notes)
    report(end+1) = "  (none)";
else
    for i = 1:numel(run_notes)
        report(end+1) = "  " + run_notes(i);
    end
end

fprintf('%s\n', report) %the one place this script prints a summary


%% Save results

%results only (no spkcount / rate / lfp_s1: they're big, re-run the loading sections to get them back),
%plus every Setup value in `settings` so animals can be pooled and checked for matching settings later
disp('Saving results')

%create results folder
if isfolder([save_directory '\StochRes\' animal]) == false
    mkdir([save_directory '\StochRes\' animal])
end

%every Setup value (the ones that exist) in one struct
setnames = {'animal','region','grpnames','condis_grp','tr_conditional_use','spike_dir','lfp_dir','probeinfo_dir', ...
    'stim_onset_30k','stim_onset_lfp','fs_mua','win_ms','base_ms','bin_ms','smooth_ms', ...
    'resp_ms','pre_resp_ms','pre_end_ms','fake_ms','tmpl_trim','tmpl_iter','chan_z_min','resp_null_prct','keep_qc','qc_plotsave', ...
    'psd_ms','psd_nw','psd_hz','line_out','line_hz', ...
    'do_align','align_ms','max_lag_ms','n_perm', ...
    'do_S1','s1_region','lfp_load_ms','s1_shank_label','s1_ref_row','tf_win_ms','tf_step_ms','tf_ms','tf_nw', ...
    'shape_block_ms','n_shuf','shuf_p','clust_p','load_ms'};
settings = struct();
for i = 1:numel(setnames)
    if exist(setnames{i},'var')
        settings.(setnames{i}) = eval(setnames{i}); %copy the variable of that name
    end
end
settings.run_date = string(datetime('now','Format','yyyy-MM-dd HH:mm'));

%what to save (only the ones this run made, so switched-off parts don't stop the save)
savelist = {'settings','report','run_notes','condis','condgrp','filesused','ntrials','trinfo', ...
    'prb','pmap','chanlist','chrow','chshank','t_bin', ...
    'npair','pairs','pairnames', ...
    'respchans','resp','respfrac','rate_mean','rate_sem', ...
    'psd_f','psd_tot','psdm_tot', ...
    'align_methods','align_use_ms','lags','lags_null','lagpk','lagsd','lagsd_null', ...
    'al_names','tc_trial','tc_cond','tc_cond_null','tc_p', ...
    'prb_s1','pmap_s1','chanlist_s1','s1_refch','tf_t','tf_f','v1tr','ntr_s1', ...
    'r_shape_ch','shape_zpix','shape_ppix','shape_sig','shape_clustthr'};
%Quality Control details: how the response trials and responsive channels were found (keep_qc)
if keep_qc == 1
    savelist = [savelist, {'chan_medz','chan_medz_cond','respchans_cond','respscore','respnull','respthr','resp_tmpl','tmpl_rounds', ...
        'align_ref','align_ref_r'}];
end
haveit = false(1,numel(savelist));
for i = 1:numel(savelist)
    haveit(i) = exist(savelist{i},'var') == 1;
end
savelist = savelist(haveit);

fname = sprintf([animal '-StochRes_MUA']);
save([save_directory '\StochRes\' animal '\' fname], savelist{:}, '-v7.3');
disp(['saved: ' save_directory '\StochRes\' animal '\' fname '.mat'])
