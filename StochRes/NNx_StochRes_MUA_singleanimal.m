% NNx_StochRes_MUA_singleanimal.m
%
% Description: a look for stochastic resonance (SR) in V1 multi-unit activity (MUA), one animal
%   at a time, narrowed to the measures that point most directly at SR. Two condition groups (e.g. light
%   alone vs light + whisker, where the whisker input could act as the "noise"), each with one or more
%   conditions, paired by position (L_8 with LW_8). A group's conditions are never pooled.
%     1) MUA (0/1 spike onsets, 30 kHz) -> 1 ms spike counts -> Gaussian-smoothed spike rate per trial
%     2) response trials (yes/no) per condition (resp_method 1 or 2 + the spike-count rule),
%        with the optional Quality Control pictures. Only response trials are used from here on
%     3) trial alignment: one lag per trial (shared by all channels, +/- max_lag_ms) from cross-correlating
%        with the most typical response trial, plus the same on the pre-stim baseline as a null. Then two
%        probe-layout figures of the condition-mean spike rate, each pair on the same axes (group 1 black,
%        group 2 red): unaligned, and aligned
%     4) trial-to-trial correlation within each condition, unaligned vs aligned, with stats:
%        - does aligning raise it more than it does on the pre-stim null (the gain test)?
%        - each condition vs every other condition, within the unaligned and within the aligned trials
%     5) (do_S1 = 1) V1 spike-count spectra next to S1 LFP spectra (spec_ms, unshifted trials), coup_groups
%        conditions, as a grid of responsive V1 channel x every S1 channel small plots in the frequency domain
%     6) (do_S1 = 1) V1 spikes x S1 LFP coupling for every responsive V1 channel x S1 channel pair, in the bands
%        (alpha 8-15 Hz, low gamma 25-45 Hz), unshifted trials (the alignment is only for steps 3-4): spike-field coherence, amplitude-envelope correlation and
%        (do_mi = 1) Gaussian-copula mutual information, each against trial shuffles with max-over-pairs correction
%     7) post-run report (printed once at the end and saved) and the results saved with all settings
%   (tried earlier and dropped: the sliding-window
%   time-frequency shape maps, the probe-layout spectra (psd_tot) and the alignment-method picture)
%
% Inputs:
%   User-set config (Setup section): animal, region, group names + condition lists, windows (ms),
%     smoothing, alignment, spectrum, band and coupling settings, toggles
%   Spiking files picked per condition with uipickfiles,
%   <spike_dir>\<animal>\<animal>-<stim>-spiking_results.mat (saved -v7.3, read with matfile):
%     stim_spike_stimchunks (double 0/1, trials x samples x channels) - 1 at the first sample of each
%       threshold crossing, 30 kHz, 450000 samples per trial (-5 s to +10 s)
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1)
%     tr_remove (0/1, 1 x nTrials) - mask over tr_keep; only tr_keep(tr_remove == 0) are used
%     tr_remove_conditional (table, optional) - Name / TrialIdx named trial sets (NNx_TrialRemove_Conditional.m)
%   <lfp_dir>\<animal>\<animal>-<stim>-LFP.mat (do_S1 = 1, coup_groups conditions, S1 channels only):
%     stim_lfp_stimchunks (double, trials x samples x channels) - peri-stimulus LFP, 1 kHz, uV, plus the
%     same tr_keep / tr_remove
%   <probeinfo_dir>\<animal>-ProbeInfo.mat, ProbeInfo fields:
%     .Areas (cell of char, 1 x probes) - region per probe; region / s1_region pick the probes
%     .ProbeMaps (cell, 1 x probes) - [rows x shanks] raw channel IDs, 0 = no site, row 1 = top
%     .Ch_Remove (cell, 1 x probes) - raw channel IDs to leave out
%     .ShankLabels (cell, 1 x probes) - struct per probe; .Light (1 x shanks) picks the channel shown in the
%       every-trial Quality Control picture (plot_respcheck = 1 only)
%   INFERRED DATA CONTRACTS (not stated in the files, same as the other scripts that read them):
%     - stim onset is sample 150001 of stim_spike_stimchunks (30 kHz, pre = 5 s) and sample 5001 of
%       stim_lfp_stimchunks (1 kHz)
%     - channel index (dim 3) = raw channel ID (true when ChanIds = 1:Chans), both probes in one file
%     - the spiking and LFP files of a stim share tr_keep / tr_remove (checked when the LFP is loaded)
%     - 1 ms count bins and 1 kHz LFP samples line up one to one (bin 0 ms = LFP sample 0 ms), so steps 5-6
%       need bin_ms = 1 (checked)
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
%     condition, of the group 2 condition]; pairnames (cell, 1 x pairs) - e.g. 'L_8 vs LW_8'
%   t_bin (double, 1 x bins) - ms (vs stim onset) where each count bin starts; load_ms (1 x 2) - ms loaded
%   spkcount (cell, 1 x conditions) - (single, trials x bins x channels) spike counts per bin_ms bin,
%     good trials only, files of a condition pooled, channels in chanlist order (not saved)
%   trinfo (cell, 1 x conditions) - (double, trials x 2) file number within the condition, trial (tr_keep value)
%   ntrials (double, 1 x conditions) - good trials per condition
%   rate (cell, 1 x conditions) - (single, trials x bins x channels) spikes/s, spkcount smoothed with a
%     Gaussian kernel smooth_ms wide (total width) (not saved)
%   filesused (string, files x 4) - condition, file name, good trials, trials left out by tr_conditional_use
%   resp (cell, 1 x conditions) - (logical, trials x 1) true = response trial; respchans (logical, channels x 1)
%     - responsive in at least one condition; respfrac (double, 1 x conditions) - fraction of response trials
%   (the response-detection details chan_medz_cond, respchans_cond, chan_medz, respscore, respnull, respthr,
%   resp_tmpl, tmpl_rounds, chpeakz, resp_shankch, resp_detect, respcount, resp_lowcount, resp_mincount are made
%   always, but only saved with keep_qc = 1; see the comments in the response-trial section for their layout)
%   rate_mean, rate_sem (cell, 1 x conditions) - (single, bins x channels) spikes/s, unaligned mean and SEM of
%     rate over the response trials
%   align_methods (cell, 1 x 1) - 'typical'; align_ref, align_ref_r (1 x conditions) - reference trial (row of
%     trinfo{con}) and its r with the others' mean (saved with keep_qc = 1)
%   lags, lags_null (cell, 1 x 1){1 x conditions} - (double, response trials x 1) ms, one shift per trial (real
%     window / pre-stim null window), |lag| <= max_lag_ms. Aligned trial(t) = trial(t + lag)
%   lagpk (cell, 1 x 1){1 x conditions} - (double, response trials x 1) peak cross-correlation r
%   lagsd, lagsd_null (double, 1 x conditions) - ms, SD of the lags (trial-to-trial jitter)
%   align_use_ms (double, 1 x 2) - ms, the alignment window actually used (align_ms trimmed to the loaded data)
%   rate_mean_al, rate_sem_al (cell, 1 x conditions) - (single, plot bins x channels) spikes/s, aligned mean and
%     SEM over the response trials, in plot_ms; t_plot (double, 1 x plot bins) - ms
%   al_names (cell, 1 x 2) - 'unaligned', 'typical'
%   tc_trial, tc_trial_null (cell, 1 x versions){1 x conditions} - (double, response trials x channels) each
%     trial's mean r (smoothed rate time courses, Fisher z mean, back to r) with the other trials of its
%     condition, in align_use_ms / in the same-length pre-stim null window (aligned with lags_null)
%   tc_cond, tc_cond_null (cell, 1 x versions) - (double, channels x conditions) mean pairwise trial r within
%     each condition, real / null window
%   tc_gain_p, tc_gain_q (double, channels x conditions) - paired signrank over trials: aligning's gain in r
%     (aligned - unaligned) in the real window vs the same gain in the null window; q = Benjamini-Hochberg over
%     the channels of a condition
%   tc_gain (double, channels x conditions) - median over trials of (real gain - null gain)
%   tc_alvsun_p (double, channels x conditions) - paired signrank, aligned vs unaligned r (real window only)
%   tc_gainch (double, channels x conditions) - gain in r from aligning: tc_cond{typical} - tc_cond{unaligned}
%   tc_gainpair_p, tc_gainpair_q, tc_gainpair_d (double, 1 x pairs) - group 2 vs group 1 condition of each pair:
%     paired signrank over the responsive channels of tc_gainch, Benjamini-Hochberg q over the pairs, median difference
%   cmppairs (double, comparisons x 2) - every pair of conditions [con A, con B]; cmpnames (cell, 1 x comparisons)
%   tc_cmp_p, tc_cmp_q (cell, 1 x versions) - (double, channels x comparisons) label-shuffle permutation p of
%     con B - con A (only those two conditions' trials shuffled, n_perm times); q = Benjamini-Hochberg over all
%     channels x comparisons of that version
%   tc_p (cell, 1 x versions) - (double, channels x pairs) the columns of tc_cmp_p for the group 1 vs group 2 pairs
%   (do_S1 = 1)
%   prb_s1, pmap_s1, chanlist_s1, chrow_s1, chshank_s1 - the S1 probe, as prb / pmap / chanlist / chrow / chshank
%   lfp_s1 (cell, 1 x conditions) - (single, response trials x samples x S1 channels) uV, S1 LFP in
%     lfp_load_ms (1 kHz), same trials/order as the MUA, coup_groups conditions only (empty otherwise; not saved)
%   t_lfp (double, 1 x samples) - ms of each LFP sample
%   coupcons (double, 1 x conditions used) - the conditions that got steps 5-6 (coup_groups, >= 3 response trials)
%   spec_use_ms (double, 1 x 2) - ms, spec_ms trimmed to the loaded data
%   spec_f (double, freqs x 1) - Hz, psd_hz range (multitaper, nfft steps)
%   v1psd (cell, 1 x conditions) - (double, freqs x channels) V1 spike-count spectrum, mean over the
%     response trials of each trial's spectrum / its random (Poisson) firing level: 1 = no rhythm
%   s1psd (cell, 1 x conditions) - (double, freqs x S1 channels) dB, S1 LFP power, mean over the same trials
%   s1apq (cell, 1 x conditions) - (double, S1 channels x 3) knee-model 1/f fit of s1psd: b, sqrt(k), exponent
%     (log10 power = b - log10(k + f^exponent), fit over psd_hz minus line_hz, refit under the first fit, as
%     NNx_LFPCorrelation_LvsLW.m); knee frequency = k^(1/exponent)
%   s1bg, s1osc (cell, 1 x conditions) - (double, freqs x S1 channels) dB, that background and s1psd above it (0 = no rhythm)
%   rch (double, responsive V1 channels x 1) - places in chanlist of the responsive V1 channels (rows of the pair
%     results, in chanlist order)
%   coh, envr, mi (cell, 1 x conditions) - (double, responsive V1 ch x S1 ch x bands) spike-field coherence
%     (band mean), envelope correlation (Fisher-z mean over trials, back to r), mutual information (bits,
%     do_mi = 1); empty for conditions not in coupcons
%   coh_z, envr_z, mi_z (cell, 1 x conditions) - (double, same size) z against the n_shuf trial shuffles
%   coh_sig, envr_sig, mi_sig (cell, 1 x conditions) - (logical, same size) beats the (1 - pair_p) percentile of
%     the biggest shuffle z over all pairs (coherence, MI: one-sided; envelope r: two-sided, |z|)
%   coup_thr (double, 3 measures x bands x conditions) - those thresholds (coh; envr; mi)
%   run_notes (string) - things worth knowing, for the report
%
% Dependencies: NNx_BaseAnalysis*.m (spiking + LFP files, ProbeInfo); NNx_TrialRemove_Conditional(_callable).m
%   (tr_remove mask format, optional conditional sets); NNx_editProbeInfo_ChLabel_StimData.m (ShankLabels.Light,
%   only for plot_respcheck); uipickfiles.m; Signal Processing Toolbox (xcorr, dpss, butter, zp2sos,
%   filtfilt, hilbert); Statistics Toolbox (corr, prctile, signrank, tiedrank, norminv, mad);
%   Parallel Computing Toolbox optional (parfor); MATLAB R2022a+ (clim, tiledlayout 'tight', boxchart, pagemtimes)
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
% condis_grp{1} = {'L_8','L_8_TTX'}; %conditions in group 1
% condis_grp{2} = {'LW_8', 'LW_8_TTX'}; %conditions in group 2

%leave out trials flagged under these conditional trial set names? (made with NNx_TrialRemove_Conditional.m,
%e.g. {"whisker_twitch_artifact"}). {} = no, only tr_remove
tr_conditional_use = {};

%where are the spiking, LFP and ProbeInfo files?
spike_dir = 'E:\Roy\Processed Silicon Probe Data\Spiking';
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%Where should I, Matlab, put the results? (Desktop while we're testing)
save_directory = 'C:\Users\rmissa\Desktop\StochRes_test'; %results in \StochRes\<animal>, figures in \Figures_StochRes\<animal>
% save_directory = 'E:\Roy\Processed Silicon Probe Data\StochRes';

%%%Analysis details%%%
stim_onset_30k = 150001; %sample of stimulus onset (0 ms) in stim_spike_stimchunks, 30 kHz
stim_onset_lfp = 5001; %sample of stimulus onset (0 ms) in stim_lfp_stimchunks, 1 kHz
fs_mua = 30000; %Hz, MUA sampling rate
win_ms = [0 3000]; %ms post stimulus, the response window everything is measured in
base_ms = [-5000 0]; %ms, pre-stim stretch loaded (fake windows for the response trials + alignment null). The files hold 5000 ms pre-stim, no more
pad_ms = 500; %ms loaded past the last spectrum window, so the band-pass filters and the Hilbert transform have no edge right at the window's end
bin_ms = 1; %ms, spike count bin size (1 ms -> spike rate and spectra at 1 kHz, the same as the LFP)
smooth_ms = 25; %ms, total width of the Gaussian kernel for the spike-rate traces (not the SD; ~ SD x 5)

%%%Response trial details%%%
%only trials with a response are used from the trial plots on (no-response trials just muddy things).
resp_ms = [0 2500]; %ms post stimulus, the stretch the response has to be in (the full response, P1 + P2)
pre_resp_ms = 400; %ms, length of the "before" stretch whose median is each trial's baseline (method 1)
pre_end_ms = 0; %ms vs the (real or fake) stimulus, where the "before" stretch ends (method 1)
fake_ms = -4800:250:-2000; %ms, fake stimulus times in the pre-stim baseline (method 1). Needs fake_ms(1) + pre_end_ms - pre_resp_ms >= base_ms(1) and fake_ms(end) + resp_ms(2) <= 0
tmpl_trim = 20; %percent of trials left out (half each end) at every time point when averaging the template (method 1)
tmpl_iter = 2; %rounds of rebuilding the template from the trials that passed (method 1)
chan_z_min = 1; %a channel counts as responsive (in a condition) when its typical (median) trial z is at least this (method 1)
resp_null_prct = 99; %a trial is a response trial when its score beats this percentile of its condition's fake scores (method 1)

%which way to find the response trials?
resp_method = 2;
    % 1 = template matching against fake no-stim windows (above; struggled with variable responses)
    % 2 = per channel: does the rate rise far above the pre-stim baseline of ALL trials (robust: median / MAD)
    %     and stay there? Response = channels 1..resp_shank_nch of one shank all do, together
base_use_ms = [-2000 0]; %ms, pre-stim stretch of every trial the baseline comes from (method 2)
peak_win_ms = 100; %ms, the rate is averaged over this sliding window before looking for the peak (method 2)
resp_z = 3; %a channel responds in a trial when it is at least this many robust SDs (1.4826 x MAD) above the baseline median (method 2; also picks the responsive channels)
resp_shank_nch = 3; %response = the first this many KEPT channels (from the top) of any one shank all respond (method 2)
resp_min_ms = 150; %ms, each of those channels has to STAY >= resp_z for at least this long (no single blips) (method 2)
resp_together_ms = 150; %ms, and they all have to do it within this stretch of each other (method 2)

%enough spikes for the later analyses? (applied to the response trials after either method)
resp_count_ms = [350 3000]; %ms, window the spikes are counted in
resp_min_spk = 10; %spikes per responsive channel (on average, summed over the responsive channels) a response trial needs there; 0 = no rule

%%%Trial alignment details%%%
%every response trial gets one shift, found by cross-correlation with the most typical response trial
align_ms = [350 3000]; %ms post stimulus, the stretch the cross-correlation looks at (sustained/rhythmic part)
max_lag_ms = 50; %ms, biggest shift allowed for a trial (either way)
n_perm = 1000; %label shuffles for the condition vs condition trial-to-trial correlation tests
fdr_q = 0.05; %Benjamini-Hochberg false discovery rate for the trial-to-trial correlation tests

%%%S1 details%%%
%compare V1 spiking with the S1 (probe 2) LFP? only turn on for animals where probe 2 was really in S1
do_S1 = 1;
    % 0 = no, stop after the trial-to-trial correlation
    % 1 = yes, spectra grid + pair coupling
s1_region = "S1"; %ProbeInfo.Areas name of the S1 probe
lfp_load_ms = [-5000 5000]; %ms, S1 LFP stretch kept (has to hold spec_ms with room for the filters)
%which group's conditions get the S1 comparisons?
coup_groups = 2;
    % 2 = group 2 only (e.g. LW)
    % [1 2] = both groups (group 1 as the control: LW - L differences per pair get their own panels)

%%%Spectrum details%%%
spec_ms = [75 3000]; %ms post stimulus, window for the V1 / S1 spectra and the pair coupling
psd_nw = 4; %multitaper time-halfbandwidth (smooths each spectrum over ~ +/- psd_nw/window length, ~1.4 Hz for 2.9 s)
psd_hz = [3 80]; %Hz, range to look at (from 3: below that the window only holds a few cycles)
%leave out line noise from the spectra? (only needed if 50 Hz mains shows up in the data)
line_out = 0;
    % 0 = no, keep all frequencies
    % 1 = yes, leave out line_hz
line_hz = [48 52]; %Hz, 50 Hz mains band (only used when line_out = 1). One row per band, e.g. [48 52; 98 102]

%%%Pair coupling details%%%
bands = [8 15; 25 45]; %Hz, one row per band
bandnames = {'alpha','lowgamma'}; %a name per band (figure titles, file names)
env_step_ms = 25; %ms, the band envelopes / analytic signals are taken every this many ms (neighbouring samples less alike)
%mutual information as well? (Gaussian-copula MI of the band signals, phase and amplitude together; slower)
do_mi = 0;
    % 0 = no, coherence + envelope correlation only
    % 1 = yes, add MI
n_shuf = 500; %trial shuffles for the pair null (V1 of one trial with S1 of another)
pair_p = 0.05; %a pair counts when it beats the (1 - pair_p) percentile of the biggest chance z over all pairs

%%%Quality Control%%%
%want to dig into how the early steps handled the data? Everything for that is switched on/off here
%responsive channels on the probe? (each channel's best typical z over the conditions, responsive ones circled)
plot_respchans = 1;
    % 0 = no
    % 1 = yes
%every trial of each condition as a heatmap, with which got called a response + mean response vs no-response traces,
%plus each condition's response template?
plot_respcheck = 1;
    % 0 = no
    % 1 = yes (one figure per condition + one template figure)
%save the Quality Control figures? (own subfolder: <save_directory>\Figures_StochRes\<animal>\QC)
qc_plotsave = 1;
    % 0 = no, just look at them
    % 1 = yes
%keep the response detection details? (saved with the results + listed in the post-run report)
keep_qc = 0;
    % 0 = no: only the trial calls (resp), responsive channels (respchans) and response fractions are kept
    % 1 = yes: also the per-trial scores, thresholds, per-condition channel z and responsive channels, the
    %     templates and the alignment reference trial of each condition

%%%Plotting details%%%
plot_ms = [-500 3000]; %ms, time range shown in the trial and condition-mean plots
%condition means in the probe layout, unaligned and aligned? (each pair on the same axes)
plot_means = 1;
%trial-to-trial correlation pictures?
plot_corr = 1;
%V1 spike spectra x S1 LFP spectra grid? (do_S1 = 1; lots of small plots, takes a while to draw)
plot_psdgrid = 1;
%pair coupling heatmaps? (do_S1 = 1)
plot_coup = 1;
%save the figures? (as PNGs; a figure with a slider gets one PNG per panel)
plotsave = 1;
png_res = 150; %dpi of the saved PNGs (a full-screen figure -> ~2900 x 1600 px)

%delete big in-between variables as we go? (saves memory)
mr_clean = 1;

%%%Speed%%%
%run the pair coupling shuffles (one S1 channel per worker) on several CPU cores?
%Needs the Parallel Computing Toolbox; without it everything runs on one core as normal
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
% load_ms = [min([base_ms(1) plot_ms(1)]) max([win_ms(2) plot_ms(2)])]; %commented out: the spectrum window plus room for the band-pass filters (pad_ms) has to fit too
load_ms = [min([base_ms(1) plot_ms(1)]) max([win_ms(2) plot_ms(2) spec_ms(2)+pad_ms])]; %ms around stim onset
%load_ms =[min([base_ms(1) plot_ms(1)])-max_lag_ms, max([win_ms(2) plot_ms(2)])+max_lag_ms]; %commented out: padding meant reloading for every max_lag_ms change. The alignment trims its window to fit instead
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
%
%resp_method = 2 (template matching above struggled with responses that vary a lot from trial to trial):
%  1) each channel's smoothed rate, averaged over a sliding peak_win_ms window
%  2) its baseline = the pre-stim stretch base_use_ms of ALL the condition's trials together: median, and
%     spread = 1.4826 x MAD (robust, so bursts in the baseline don't inflate it; 1.4826 x MAD = SD for
%     normal data). Robust z = (rate - median) / spread
%  3) per trial and channel: the level it HOLDS for resp_min_ms, anywhere in resp_ms (any latency and any
%     shape counts, but single blips don't)
%  4) a trial is a response when channels 1..resp_shank_nch (kept, from the top) of ANY one shank all hold
%     >= resp_z, within resp_together_ms of each other: per shank and time the lowest of their held levels
%     (each taken as the best within resp_together_ms), its best time, per trial the best shank (respscore)
%  (tried first: each channel's single highest peak anywhere in resp_ms, at unrelated times; too permissive)
%  5) responsive channels (for the later sections): the condition's trial-mean trace peaks >= resp_z robust
%     SDs above its own baseline
%Each channel on its own (no pooling across channels: channels differ in latency and shape)
disp(['MUA: Finding the trials with a response (method ' num2str(resp_method) ')'])
tic

if resp_method == 1
    %do the fake windows fit? (inside the loaded stretch, and ending before the real stimulus)
    if fake_ms(1) + pre_end_ms - pre_resp_ms < load_ms(1)
        error(['The first fake window starts before the loaded data (' num2str(load_ms(1)) ' ms). Move fake_ms later or base_ms earlier'])
    end
    if fake_ms(end) + resp_ms(2) > 0
        error('The last fake response window runs past the real stimulus. Move fake_ms earlier')
    end
else
    %does the baseline stretch fit? (inside the loaded stretch, before the stimulus)
    if base_use_ms(1) < load_ms(1) || base_use_ms(2) > 0
        error(['base_use_ms has to lie between the start of the loaded data (' num2str(load_ms(1)) ' ms) and 0'])
    end
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
tmpl_rounds = zeros(1,ncond); %template rounds actually run per condition (0 = method 2, no template)

if resp_method == 2
    baseidx = find(t_bin >= base_use_ms(1) & t_bin < base_use_ms(2)); %bins of the baseline stretch
    pkw = max(1, round(peak_win_ms/bin_ms)); %bins per peak window

    %each shank's first resp_shank_nch kept channels (from the top): a response needs all of them on one shank
    nshk = size(pmap,2);
    resp_shankch = nan(nshk, resp_shank_nch); % [shanks x channels] raw IDs (NaN row = shank not used)
    shkidx = cell(1,nshk); %{shanks} their places in chanlist
    for shank = 1:nshk
        kept = pmap(:,shank); %top to bottom
        kept = kept(kept > 0 & ismember(kept, chanlist)); %kept channels only
        if numel(kept) < resp_shank_nch
            run_notes(end+1) = "Shank " + shank + ": only " + numel(kept) + " kept channel(s), not used for the response calls (resp_shank_nch = " + resp_shank_nch + ")";
            continue
        end
        resp_shankch(shank,:) = kept(1:resp_shank_nch)';
        [~,shkidx{shank}] = ismember(kept(1:resp_shank_nch), chanlist);
    end
    useshk = find(isnan(resp_shankch(:,1)) == false)'; %shanks used
    if isempty(useshk)
        error(['No shank has ' num2str(resp_shank_nch) ' kept channels, so no response can be called. Lower resp_shank_nch'])
    end
    shch = unique([shkidx{useshk}]); %the shanks' channels 1..resp_shank_nch (places in chanlist): the only ones the call needs
    nmin = max(1, round(resp_min_ms/bin_ms)); %bins a channel has to stay up
    ntog = max(1, round(resp_together_ms/bin_ms)); %bins the channels of a shank have to do it within

    chpeakz = cell(1,ncond); %{conditions} [trials x channels] each channel's peak robust z in resp_ms
    for con = 1:ncond
        ntr = ntrials(con);
        disp(['MUA: Peaks above the all-trial baseline, ' condis{con}])

        W = movmean(rate{con}, pkw, 2); % [trials x bins x channels] spikes/s, averaged over peak_win_ms

        %each channel's robust baseline over ALL trials x ALL baseline bins
        B = reshape(W(:,baseidx,:), ntr*numel(baseidx), nch); % [(trials*baseline bins) x channels]
        bmed = median(B,1); % [1 x channels]
        bsd = 1.4826 * mad(B,1,1); % [1 x channels] robust SD (median absolute deviation, scaled)
        nomad = bsd == 0;
        if any(nomad) %very sparse channel: MAD is 0, use the SD instead
            bsd(nomad) = std(B(:,nomad),0,1);
            run_notes(end+1) = string(condis{con}) + ": baseline MAD = 0 on CH" + strjoin(string(chanlist(nomad)'),", CH") + ", SD used instead";
        end
        bsd(bsd == 0) = NaN; %silent channel, drops out (and its shank can't pass)
        % Z = (W - reshape(bmed,1,1,nch)) ./ reshape(bsd,1,1,nch); %commented out: only the shanks' channels 1..3 are needed now
        Z = (W(:,:,shch) - reshape(bmed(shch),1,1,[])) ./ reshape(bsd(shch),1,1,[]); % [trials x bins x those channels] robust z, each vs ITS OWN baseline

        %commented out: one peak anywhere in resp_ms per channel, at any time, was too permissive (single blips, and
        %the three channels' peaks could fall at unrelated times)
        % chpk = reshape(max(Z(:,respidx,:),[],2), ntr, nch);
        % chpkn = reshape(max(Z(:,baseidx,:),[],2), ntr, nch);
        % chpeakz{con} = chpk;
        %tighter: each channel has to STAY up for resp_min_ms, and channels 1..3 of a shank have to do so TOGETHER
        Zs = movmin(Z, nmin, 2); %the level each channel holds for a whole resp_min_ms stretch around each time
        Zt = movmax(Zs, ntog, 2); %the best such level within resp_together_ms around each time
        chpeakz{con} = nan(ntr,nch); %each channel's sustained peak in resp_ms (only the shanks' channels 1..3)
        chpeakz{con}(:,shch) = reshape(max(Zs(:,respidx,:),[],2), ntr, numel(shch));

        %trial call: per shank, at each time the level ALL its channels 1..resp_shank_nch hold (lowest of them; a
        %silent channel = NaN = the shank can't pass), its best time in resp_ms; per trial the best shank
        shk = nan(ntr,nshk); % [trials x shanks]
        shkn = nan(ntr,nshk);
        for shank = useshk
            % shk(:,shank) = min(chpk(:,shkidx{shank}),[],2); %commented out with the per-channel peaks above
            % shkn(:,shank) = min(chpkn(:,shkidx{shank}),[],2);
            [~,loc] = ismember(shkidx{shank}, shch); %this shank's channels in Z
            together = min(Zt(:,:,loc),[],3,'includenan'); % [trials x bins] level all of them hold, close together in time
            shk(:,shank) = max(together(:,respidx),[],2);
            shkn(:,shank) = max(together(:,baseidx),[],2);
        end
        respscore{con} = max(shk,[],2,'omitnan'); % [trials x 1] level the best shank's channels 1..3 all hold for resp_min_ms, together
        respnull{con} = max(shkn,[],2,'omitnan'); % [trials x 1] the same score in the baseline stretch (QC only)
        respthr(con) = resp_z;
        resp{con} = respscore{con} >= resp_z;

        %responsive channels (later sections): the trial-mean trace's peak in resp_ms vs its own baseline, robust z
        M = reshape(mean(W,1), size(W,2), nch); % [bins x channels]
        mmed = median(M(baseidx,:),1);
        msd = 1.4826 * mad(M(baseidx,:),1,1);
        msd(msd == 0) = NaN;
        chan_medz_cond(:,con) = ((max(M(respidx,:),[],1) - mmed) ./ msd)';
        respchans_cond(:,con) = chan_medz_cond(:,con) >= resp_z;

        %the "template" for the QC picture: mean of the response trials above the baseline median
        if sum(resp{con}) >= 1
            resp_tmpl{con} = reshape(mean(rate{con}(resp{con},respidx,:),1), nw, nch) - bmed;
        else
            resp_tmpl{con} = nan(nw,nch,'single');
        end
    end
    clear W B Z Zs Zt M together shk shkn

else %resp_method == 1, template matching against fake windows
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
end

%channels responsive in any condition (what the later sections use: alignment, plot labels)
chan_medz = max(chan_medz_cond,[],2); % [channels x 1] best typical trial z over the conditions
respchans = any(respchans_cond,2); % [channels x 1]
if any(respchans) == false
    if resp_method == 2
        error(['No channel''s condition-mean trace peaks >= resp_z (' num2str(resp_z) ') robust SDs above its baseline in any condition. Try a lower resp_z'])
    else
        error(['No channel has a typical trial z >= chan_z_min (' num2str(chan_z_min) ') in any condition, so no responses can be called. Try a lower chan_z_min'])
    end
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

%enough spikes for the later analyses? A detected response trial is only kept when its raw spike count in
%resp_count_ms, summed over the responsive channels (the ones the later sections use), reaches
%resp_min_spk x the number of responsive channels. Same rule for every condition (no group rate bias)
cntidx = find(t_bin >= resp_count_ms(1) & t_bin < resp_count_ms(2)); %bins counted
resp_detect = resp; %{conditions} [trials x 1] what detection called, before the spike-count rule
respcount = cell(1,ncond); %{conditions} [trials x 1] spikes in resp_count_ms, summed over the responsive channels
resp_lowcount = zeros(1,ncond); %detected trials dropped per condition for too few spikes
resp_mincount = resp_min_spk * sum(respchans); %spikes a response trial needs
for con = 1:ncond
    respcount{con} = reshape(sum(sum(spkcount{con}(:,cntidx,respchans),2),3), ntrials(con), 1);
    if resp_min_spk > 0
        keep = respcount{con} >= resp_mincount;
        resp_lowcount(con) = sum(resp{con} & keep == false);
        resp{con} = resp{con} & keep;
    end
end

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

%per condition: every trial as a heatmap (rate of ONE V1 channel: kept channel 2 of the shank with the highest
%Light score; recording order) and the same trials as stacked lines, a strip marking which got called a response, and the
%mean of the response vs no-response trials underneath. Then one figure of every condition's final response
%template (summed over its responsive channels)
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

    %the channel shown: always the 2nd KEPT channel (from the top) of the V1 shank with the highest "Light" score
    %(ProbeInfo.ShankLabels); fixed on purpose, not a Setup option
    v1_shank_label = "Light";
    v1_ref_row = 2;
    if isfield(ProbeInfo,'ShankLabels') == false || numel(ProbeInfo.ShankLabels) < prb || isstruct(ProbeInfo.ShankLabels{prb}) == false ...
            || isfield(ProbeInfo.ShankLabels{prb}, char(v1_shank_label)) == false
        haslabels = 'none';
        if isfield(ProbeInfo,'ShankLabels') && numel(ProbeInfo.ShankLabels) >= prb && isstruct(ProbeInfo.ShankLabels{prb})
            haslabels = strjoin(fieldnames(ProbeInfo.ShankLabels{prb}),', ');
        end
        error(['No ShankLabels.' char(v1_shank_label) ' on the V1 probe. Labels there: ' haslabels])
    end
    [~,pkv] = max(ProbeInfo.ShankLabels{prb}.(char(v1_shank_label))); %shank with the highest score
    keptv = pmap(:, pkv); %that shank, top to bottom
    keptv = keptv(keptv > 0 & ismember(keptv, chanlist)); %its kept channels only, still top to bottom
    if numel(keptv) < v1_ref_row
        error(['V1 shank ' num2str(pkv) ' has only ' num2str(numel(keptv)) ' kept channels, v1_ref_row is ' num2str(v1_ref_row)])
    end
    v1_qcch = keptv(v1_ref_row); %raw channel ID
    chi = find(chanlist == v1_qcch); %its place in chanlist

    for con = 1:ncond
        %trsum = sum(rate{con}(:,pidx,respchans),3); %commented out: responsive channels are per condition now
        %trsum = sum(rate{con}(:,pidx,respchans_cond(:,con)),3); %commented out: the sum over channels hid single-channel responses, one channel now
        trsum = reshape(rate{con}(:,pidx,chi), ntrials(con), numel(pidx)); % [trials x bins shown] spikes/s, channel v1_qcch
        rtr = resp{con} == true;
        if respchans_cond(chi,con) == true
            rtxt = 'responsive';
        else
            rtxt = 'not responsive';
        end

        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        % tl = tiledlayout(3,12,'TileSpacing','compact','Padding','compact'); %commented out: a second column (stacked trials) was added
        tl = tiledlayout(3,24,'TileSpacing','compact','Padding','compact'); %left half: heatmap, strip, means; right half: every trial stacked

        %all trials as an image, in recording order
        ax1 = nexttile(tl,1,[2 11]);
        imagesc(ax1, t_plot, 1:ntrials(con), trsum)
        xline(ax1,0,'w')
        xline(ax1,resp_ms,'w--')
        cb = colorbar(ax1);
        cb.Label.String = ['spikes/s, CH' num2str(v1_qcch)];
        ylabel(ax1,'trial (recording order)')
        set(ax1,'XTickLabel',[])
        title(ax1,{[animal ' ' condis{con} ': ' num2str(sum(rtr)) ' of ' num2str(ntrials(con)) ' trials called a response (call made on the responsive channels)'], ...
            ['showing CH' num2str(v1_qcch) ': kept channel ' num2str(v1_ref_row) ' of the peak ' char(v1_shank_label) ' shank (shank ' num2str(pkv) '), ' rtxt ' in this condition']},'Interpreter','none')

        %strip beside it: green = response, grey = no response, orange = detected but too few spikes (resp_min_spk)
        ax2 = nexttile(tl,12,[2 1]);
        lowcnt = resp_detect{con} == true & rtr == false; %detected, dropped by the spike-count rule
        stripval = double(rtr);
        stripval(lowcnt) = 2;
        % imagesc(ax2, 1, 1:ntrials(con), double(rtr)) %commented out: a third colour for the spike-count rule
        imagesc(ax2, 1, 1:ntrials(con), stripval)
        colormap(ax2,[0.8 0.8 0.8; 0 0.6 0; 1 0.55 0])
        clim(ax2,[0 2])
        set(ax2,'XTick',[],'YTick',[])
        title(ax2,'resp?')

        %mean trace of the response vs no-response trials
        % ax3 = nexttile(tl,25,[1 11]); %commented out: the grid is 24 wide now
        ax3 = nexttile(tl,49,[1 11]); %row 3, left half
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

        %right half: every trial as its own line, one below the other (trial 1 on top), green = response, grey = no response
        ax4 = nexttile(tl,13,[3 12]);
        off = prctile(trsum(:),95); %spikes/s each trial is moved down: most of a trial's trace stays clear of the next one
        if isnan(off) || off <= 0
            off = 1;
        end
        hold(ax4,'on')
        for tr = 1:ntrials(con)
            if rtr(tr) == true
                col = [0 0.6 0]; %response
            elseif lowcnt(tr) == true
                col = [1 0.55 0]; %detected, too few spikes
            else
                col = [.6 .6 .6]; %no response
            end
            plot(ax4, t_plot, trsum(tr,:) - (tr-1)*off, 'Color',col, 'LineWidth',0.75)
        end
        xline(ax4,0)
        xline(ax4,resp_ms,'--')
        hold(ax4,'off')
        xlim(ax4,plot_ms)
        ylim(ax4,[-(ntrials(con)-1)*off - off, max(trsum(1,:)) + off])
        set(ax4,'YTick',fliplr(-(0:ntrials(con)-1)*off),'YTickLabel',fliplr(trinfo{con}(:,2)')) %each line labelled with its trial number
        xlabel(ax4,'ms from stim onset')
        ylabel(ax4,'trial (recording order)')
        title(ax4,['every trial, CH' num2str(v1_qcch) ' (each moved down ' num2str(off,'%.0f') ' spikes/s): green = response, orange = detected but < ' num2str(resp_mincount) ' spikes in ' num2str(resp_count_ms(1)) '-' num2str(resp_count_ms(2)) ' ms, grey = no response'])

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



%Now beginning the S1 loading (needed by the spectra grid and the pair coupling at the end)

%% Data Extraction: Load the S1 LFP of the same recordings

%the LFP files sit next to the spiking files (...\LFP\<animal>\<animal>-<stim>-LFP.mat), so no picking.
%Only the S1 probe's kept channels, only the coup_groups conditions, same response trials in the same order as the MUA
%(an earlier version also picked a single S1 reference channel and shank for its maps; not needed any more)
if do_S1 == 1
    disp('S1: Loading the S1 LFP of the same recordings')
    tic

    if bin_ms ~= 1
        error('The V1 - S1 comparisons need bin_ms = 1 (spike bins = LFP samples, 1 kHz)')
    end
    lfpsamp = stim_onset_lfp + (lfp_load_ms(1):lfp_load_ms(2)-1); %LFP samples to keep (1 kHz)
    t_lfp = lfp_load_ms(1):lfp_load_ms(2)-1; %ms of each kept sample

    %S1 probe channels, and where they sit on the probe
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
    nch_s1 = numel(chanlist_s1);
    chrow_s1 = zeros(nch_s1,1); %ProbeMaps row of each kept S1 channel (1 = top)
    chshank_s1 = zeros(nch_s1,1); %ProbeMaps column (shank)
    for c1 = 1:nch_s1
        [chrow_s1(c1), chshank_s1(c1)] = find(pmap_s1 == chanlist_s1(c1));
    end

    lfp_s1 = cell(1,ncond); %{conditions} [response trials x samples x S1 channels] uV, single (empty = not a coup_groups condition)
    for con = find(ismember(condgrp, coup_groups))
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



%% Data Analysis: Condition-mean spike rate of the aligned trials

%each response trial read at t + its lag (the aligned version of rate_mean / rate_sem), in plot_ms. A shift that
%would run off the loaded data is clipped at the edge (it can't with the default windows: load_ms has max_lag_ms spare)
pidx = find(t_bin >= plot_ms(1) & t_bin < plot_ms(2)); %bins shown
t_plot = t_bin(pidx); %ms
nplot = numel(pidx);
rate_mean_al = cell(1,ncond); %{conditions} [plot bins x channels] spikes/s, aligned response trials
rate_sem_al = cell(1,ncond); %{conditions} [plot bins x channels] spikes/s
for con = 1:ncond
    rtrs = find(resp{con} == true); %response trials
    nresp = numel(rtrs);
    lb = round(lags{1}{con} / bin_ms); %bins, each trial's shift
    r_al = zeros(nresp,nplot,nch,'single'); % [response trials x plot bins x channels]
    for tr = 1:nresp
        idx = min(max(pidx + lb(tr), 1), nbins); %this trial read at t + its lag
        r_al(tr,:,:) = rate{con}(rtrs(tr), idx, :);
    end
    rate_mean_al{con} = reshape(mean(r_al,1), nplot, nch);
    rate_sem_al{con} = reshape(std(r_al,0,1), nplot, nch) / sqrt(nresp);
end
clear r_al


%% Data Visualization: Condition means in the probe layout, unaligned and aligned

%two figures, the V1 probe as it sits (one small plot per site, ProbeMaps layout): unaligned, then aligned.
%Each pair's two conditions on the same axes: group 1 condition black, group 2 condition red, mean +/- SEM
%(shaded) of the response trials. Slider over the pairs. Every channel keeps the same y range in both figures,
%so you can flip between them and see what the alignment did. Bold CH label = responsive channel
if plot_means == 1
    disp('MUA: Making condition-mean pictures (unaligned and aligned)')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    [nrow,nshank] = size(pmap);
    paircol = [0 0 0; 0.85 0 0]; %group 1 condition black, group 2 condition red

    %y range per pair and channel: the biggest mean + SEM of either condition, unaligned or aligned
    ytop = zeros(npair,nch); % [pairs x channels] spikes/s
    for k = 1:npair
        for j = 1:2
            con = pairs(k,j);
            top1 = max(rate_mean{con}(pidx,:) + rate_sem{con}(pidx,:),[],1); % [1 x channels] unaligned
            top2 = max(rate_mean_al{con} + rate_sem_al{con},[],1); %aligned
            ytop(k,:) = max([ytop(k,:); top1; top2],[],1,'omitnan');
        end
    end
    ytop(ytop <= 0 | isnan(ytop)) = 1;
    ytop = ytop * 1.05;

    for fv = 1:2 %1 = unaligned, 2 = aligned
        if fv == 1
            vtxt = 'unaligned';
        else
            vtxt = ['aligned (each trial shifted by its lag, typical-trial template, ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms, max ' num2str(max_lag_ms) ' ms)'];
        end
        meanfig = figure();
        set(meanfig, 'Position', get(0, 'Screensize'));
        for k = 1:npair
            %this pair's panel (only the first one starts visible); Tag + UserData let the slider find it
            pan = uipanel(meanfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
                'Tag','condpanel', 'UserData',k, 'Visible', k == 1);
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
                    hold(ax,'on')
                    for j = 1:2 %group 1 condition, then group 2 condition
                        con = pairs(k,j);
                        if fv == 1
                            m = rate_mean{con}(pidx,ch)'; % [1 x plot bins]
                            s = rate_sem{con}(pidx,ch)';
                        else
                            m = rate_mean_al{con}(:,ch)';
                            s = rate_sem_al{con}(:,ch)';
                        end
                        if all(isnan(m)) %no response trials
                            continue
                        end
                        s(isnan(s)) = 0; %one trial: no SEM
                        fill(ax, [t_plot fliplr(t_plot)], [m-s fliplr(m+s)], paircol(j,:), 'FaceAlpha',0.2, 'EdgeColor','none') %SEM band
                        plot(ax, t_plot, m, 'Color',paircol(j,:), 'LineWidth',1.2)
                    end
                    xline(ax,0,'b')
                    hold(ax,'off')
                    xlim(ax,plot_ms)
                    ylim(ax,[0 ytop(k,ch)])
                    if respchans(ch) == true
                        text(ax,0.02,0.9,['CH' num2str(chid)],'Units','normalized','FontWeight','bold')
                    else
                        text(ax,0.02,0.9,['CH' num2str(chid)],'Units','normalized','Color',[.45 .45 .45])
                    end
                    if row < nrow
                        set(ax,'XTickLabel',[])
                    end
                end
            end
            ca = pairs(k,1);
            cb = pairs(k,2);
            xlabel(tl,'ms from stim onset')
            ylabel(tl,'spikes/s (each channel its own range, same in both figures)')
            title(tl,{[animal ' ' ProbeInfo.Areas{prb} '-Probe condition means, ' vtxt], ...
                [condis{ca} ' (black, ' num2str(sum(resp{ca})) ' response trials) vs ' condis{cb} ' (red, ' num2str(sum(resp{cb})) '): mean +/- SEM; bold = responsive channel; pair ' num2str(k) ' of ' num2str(npair)]},'Interpreter','none')
        end

        %the slider (finds the panels by their Tag, so it still works in the saved figure)
        if npair > 1
            uicontrol(meanfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
                'String','pair:', 'HorizontalAlignment','right')
            uicontrol(meanfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
                'Min',1, 'Max',npair, 'Value',1, 'SliderStep',[1 1]/(npair-1), ...
                'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
                num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
        end
        drawnow %show it now; the script doesn't wait, it carries on

        if plotsave == 1
            if fv == 1
                ftag = 'means_unaligned';
            else
                ftag = 'means_aligned';
            end
            fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' ftag]);
            %one PNG per panel (= pair): show only panel k, save its layout (not the slider), then back to panel 1
            ptags = strrep(pairnames,' ','_'); %file name tag of each panel, e.g. L_8_vs_LW_8
            pans = findobj(meanfig,'Tag','condpanel');
            for k = 1:numel(ptags)
                set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
                pk = findobj(pans,'flat','UserData',k);
                tlk = findobj(pk,'Type','tiledlayout');
                exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' ptags{k} '.png'],'Resolution',png_res)
            end
            set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
            set(findobj(meanfig,'Style','slider'),'Value',1);
        end
    end
end




%Now beginning the trial-to-trial correlation section



%% Data Analysis: Trial-to-trial correlation, unaligned vs aligned

%How alike are a condition's response trials? Every pair of response trials, their smoothed rate time courses in
%the alignment window (align_use_ms) correlated (r), Fisher z, averaged within each condition, back to r
%(as NNx_LFPCorrelation_LvsLW.m; a group's conditions are never pooled). Done on the unaligned trials and on the
%aligned ones (typical-trial template).
%Null: the same numbers on the same-length pre-stim stretch, aligned with the null lags: lining trials up always
%makes them look a bit more alike, even pure noise, so the gain from aligning is tested against this null.
%Condition vs condition: every pair of conditions (not just L vs LW), the labels of those two conditions' trials
%shuffled n_perm times, p = how often a shuffle gives a difference at least as big as the real one
disp('MUA: Trial consistency, unaligned vs aligned')
tic

al_names = {'unaligned','typical'}; %the versions compared
al_meth = [0 find(strcmp(align_methods,'typical'))]; %which alignment method's lags each version uses (0 = no shift)
nver = numel(al_names);
wbr = winbins{1}; %bins of the alignment window (from the alignment section)
wbn = winbins{2}; %bins of its pre-stim copy (null)
nwin_al = numel(wbr);

%every pair of conditions, for the condition vs condition tests
cmppairs = nchoosek(1:ncond,2); % [comparisons x 2] condition A, condition B
ncmp = size(cmppairs,1);
cmpnames = cell(1,ncmp); %{comparisons} e.g. 'LW_8 vs L_8' (B vs A: the difference is B - A)
for i = 1:ncmp
    cmpnames{i} = [condis{cmppairs(i,2)} ' vs ' condis{cmppairs(i,1)}];
end

tc_trial = cell(1,nver); %{versions}{conditions} [response trials x channels] each trial's mean r with the other trials of its condition
tc_trial_null = cell(1,nver); %{versions}{conditions} [response trials x channels] the same in the pre-stim null window
tc_cond = cell(1,nver); %{versions} [channels x conditions] mean pairwise trial r within each condition
tc_cond_null = cell(1,nver); %{versions} [channels x conditions] the same on the pre-stim window
tc_cmp_p = cell(1,nver); %{versions} [channels x comparisons] permutation p, condition B vs condition A

for v = 1:nver
    disp(['MUA: ' al_names{v}])
    tc_trial{v} = cell(1,ncond);
    tc_trial_null{v} = cell(1,ncond);

    %line every condition's response trials up (or not) and stack them
    allr = []; % [all response trials x window bins x channels] smoothed rate, real window
    allrn = []; % [all response trials x window bins x channels] smoothed rate, pre-stim window
    condlab = []; %condition of each stacked trial
    for con = 1:ncond
        rtrs = find(resp{con} == true); %response trials
        nresp = numel(rtrs);
        tc_trial{v}{con} = nan(nresp,nch);
        tc_trial_null{v}{con} = nan(nresp,nch);
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
    ntot = numel(condlab);

    %trial-to-trial correlation, per channel
    tc_cond{v} = nan(nch,ncond);
    tc_cond_null{v} = nan(nch,ncond);
    tc_cmp_p{v} = nan(nch,ncmp);
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

            %each trial with the other trials of its condition, and per condition (never pooled)
            zc = nan(1,ncond);
            for con = 1:ncond
                idx = condlab == con;
                if sum(idx) < 2
                    continue
                end
                if w == 1
                    tc_trial{v}{con}(:,ch) = tanh(mean(Z(idx,idx),2,'omitnan'));
                else
                    tc_trial_null{v}{con}(:,ch) = tanh(mean(Z(idx,idx),2,'omitnan'));
                end
                zc(con) = mean(Z(idx,idx),'all','omitnan');
            end
            if w == 1
                tc_cond{v}(ch,:) = tanh(zc);
            else
                tc_cond_null{v}(ch,:) = tanh(zc);
                continue %the permutation tests are on the real window only
            end

            %condition B vs condition A for every pair of conditions, shuffling only those two conditions' trial labels
            for i = 1:ncmp
                ca = cmppairs(i,1);
                cb = cmppairs(i,2);
                if isnan(zc(ca)) || isnan(zc(cb))
                    continue
                end
                pti = find(condlab == ca | condlab == cb); %the two conditions' trials in the stack
                Zp = Z(pti,pti);
                plab = condlab(pti); %their condition
                obs = zc(cb) - zc(ca);
                %the mean of Zp over a set of trials = (sum over the set's pairs) / (number of those pairs that
                %aren't NaN), and for every shuffle at once that is one matrix product
                ntp = numel(plab);
                IB = zeros(ntp,n_perm); % [their trials x shuffles] 1 = labelled as condition B in that shuffle
                for p = 1:n_perm
                    shlab = plab(randperm(ntp)); %shuffled labels within the two conditions
                    IB(:,p) = shlab == cb;
                end
                IA = 1 - IB; %labelled as condition A
                Wp = double(isnan(Zp) == false); %1 = a pair that counts
                Z0 = Zp;
                Z0(Wp == 0) = 0;
                nulldiff = (sum(IB .* (Z0*IB),1) ./ sum(IB .* (Wp*IB),1) - sum(IA .* (Z0*IA),1) ./ sum(IA .* (Wp*IA),1))'; % [shuffles x 1]
                tc_cmp_p{v}(ch,i) = (sum(abs(nulldiff) >= abs(obs)) + 1) / (n_perm + 1);
            end
        end
    end
end
clear allr allrn r_al r_aln X C Z Zp Z0 Wp IA IB

%the group 1 vs group 2 pairs (tc_p): their columns of tc_cmp_p
tc_p = cell(1,nver); %{versions} [channels x pairs]
for v = 1:nver
    tc_p{v} = nan(nch,npair);
    for k = 1:npair
        i = find(cmppairs(:,1) == pairs(k,1) & cmppairs(:,2) == pairs(k,2)); %(group 1's conditions come first in condis, so A = group 1)
        tc_p{v}(:,k) = tc_cmp_p{v}(:,i);
    end
end

%Benjamini-Hochberg q of the condition vs condition tests, over all channels x comparisons of a version
tc_cmp_q = cell(1,nver); %{versions} [channels x comparisons]
for v = 1:nver
    pv = tc_cmp_p{v}(:); %every test of this version
    q = nan(size(pv));
    ok = find(isnan(pv) == false);
    [ps,ord] = sort(pv(ok)); %smallest p first
    nt = numel(ps);
    qs = ps .* nt ./ (1:nt)'; %p x tests / rank
    qs = min(1, flipud(cummin(flipud(qs)))); %can't be smaller than the q of a bigger p
    q(ok(ord)) = qs;
    tc_cmp_q{v} = reshape(q, nch, ncmp);
end

%does aligning raise the consistency more than it does on the null? Per trial: gain = aligned r - unaligned r in
%the real window, against the same gain in the pre-stim window (the same trials, so a paired signrank over trials).
%(each trial's r is with the OTHER trials, so trials aren't fully independent: read the p values as a guide)
vun = find(strcmp(al_names,'unaligned'));
vstr = find(strcmp(al_names,'typical'));
tc_gain = nan(nch,ncond); % [channels x conditions] median over trials of (real gain - null gain)
tc_gain_p = nan(nch,ncond); % [channels x conditions] signrank p, real gain vs null gain
tc_gain_q = nan(nch,ncond); % [channels x conditions] Benjamini-Hochberg over the channels of a condition
tc_alvsun_p = nan(nch,ncond); % [channels x conditions] signrank p, aligned vs unaligned r (real window)
for con = 1:ncond
    if size(tc_trial{vun}{con},1) < 3 %too few trials to test
        continue
    end
    for ch = 1:nch
        greal = tc_trial{vstr}{con}(:,ch) - tc_trial{vun}{con}(:,ch); % [response trials x 1] gain, real window
        gnull = tc_trial_null{vstr}{con}(:,ch) - tc_trial_null{vun}{con}(:,ch); %gain, null window
        ok = isnan(greal) == false & isnan(gnull) == false;
        if sum(ok) < 3
            continue
        end
        tc_gain(ch,con) = median(greal(ok) - gnull(ok));
        tc_gain_p(ch,con) = signrank(greal(ok), gnull(ok));
        tc_alvsun_p(ch,con) = signrank(tc_trial{vstr}{con}(ok,ch), tc_trial{vun}{con}(ok,ch));
    end

    %Benjamini-Hochberg over this condition's channels
    pv = tc_gain_p(:,con);
    ok = find(isnan(pv) == false);
    [ps,ord] = sort(pv(ok));
    nt = numel(ps);
    qs = ps .* nt ./ (1:nt)';
    qs = min(1, flipud(cummin(flipud(qs))));
    tc_gain_q(ok(ord),con) = qs;
end

%gain from aligning per channel (plain: aligned - unaligned mean pairwise r in the response window), and L vs LW in each
%pair: paired signrank over the responsive channels (the same channels in both conditions) of gain(LW) - gain(L).
%(both conditions carry the same small rise any alignment gives, so it cancels in the comparison. Channels of one
%animal aren't fully independent: read it as a within-animal check)
tc_gainch = tc_cond{vstr} - tc_cond{vun}; % [channels x conditions] gain in r from aligning
tc_gainpair_p = nan(1,npair); % [1 x pairs] signrank p, group 2 vs group 1 condition of each pair
tc_gainpair_d = nan(1,npair); % [1 x pairs] median over responsive channels of gain(group 2) - gain(group 1)
for k = 1:npair
    ga = tc_gainch(respchans,pairs(k,1)); %group 1 condition (e.g. L_8)
    gb = tc_gainch(respchans,pairs(k,2)); %group 2 condition (e.g. LW_8)
    ok = isnan(ga) == false & isnan(gb) == false;
    if sum(ok) < 3 %too few channels to test
        continue
    end
    tc_gainpair_d(k) = median(gb(ok) - ga(ok));
    tc_gainpair_p(k) = signrank(gb(ok), ga(ok));
end
%Benjamini-Hochberg over the pairs
tc_gainpair_q = nan(1,npair);
ok = find(isnan(tc_gainpair_p) == false);
[ps,ord] = sort(tc_gainpair_p(ok));
nt = numel(ps);
qs = ps .* nt ./ (1:nt); %p x tests / rank
qs = min(1, fliplr(cummin(fliplr(qs)))); %can't be smaller than the q of a bigger p
tc_gainpair_q(ok(ord)) = qs;
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);




%% Data Visualization: Trial-to-trial correlation results

%picture 1: the gain in trial-to-trial r from aligning (aligned - unaligned), aligned conditions only, one box per
%           condition, one dot per responsive channel; group 2 vs group 1 per pair, paired signrank over the channels
%           (was probe maps per condition with a slider)
%(picture 2, the condition x condition matrix, was dropped; the comparisons are in the report)
%picture 3: each responsive channel's consistency (mean pairwise trial r), unaligned vs aligned, one dot per channel,
%           real window in the condition's colour, pre-stim null in grey, conditions side by side in their pairs
if plot_corr == 1
    disp('MUA: Making trial-to-trial correlation pictures')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    [nrow,nshank] = size(pmap);
    [rrow,rshank] = find(pmap > 0 & ismember(pmap, chanlist) == false); %removed channels
    nc = 128; %blue - white - red colour map for differences (0 = white)
    divmap = [linspace(0,1,nc)' linspace(0.3,1,nc)' ones(nc,1); ones(nc,1) linspace(1,0.2,nc)' linspace(1,0.2,nc)'];

    %colour of each condition: group 1 greys/black, group 2 reds, subgroups dark to light
    grpdark = [0 0 0; 0.6 0.1 0.1];
    grplight = [0.55 0.55 0.55; 1 0.55 0.45];
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

    %commented out (Roy: box plot instead of the gain maps; the per-trial gain-vs-null test, tc_gain_q, is still in the report)
    %{
    %picture 1: the gain maps, one panel per condition
    gainvals = {}; %{conditions}{3 maps} values (channels x 1)
    for con = 1:ncond
        gainvals{con} = {tc_cond{vstr}(:,con) - tc_cond{vun}(:,con), tc_cond_null{vstr}(:,con) - tc_cond_null{vun}(:,con), tc_gain(:,con)};
    end
    allg = [];
    for con = 1:ncond
        allg = [allg; cell2mat(gainvals{con}')];
    end
    dmax = max(abs(allg),[],'omitnan'); %one symmetric colour range for every map and condition
    if isempty(dmax) || isnan(dmax) || dmax == 0
        dmax = 0.01;
    end
    gtitles = {'gain from aligning: aligned - unaligned r', 'same gain on the pre-stim null', ['per-trial gain beyond the null, median (square = q < ' num2str(fdr_q) ')']};

    gainfig = figure();
    set(gainfig, 'Position', get(0, 'Screensize'));
    for con = 1:ncond
        pan = uipanel(gainfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
            'Tag','condpanel', 'UserData',con, 'Visible', con == 1);
        tl = tiledlayout(pan,1,3,'TileSpacing','compact','Padding','compact');
        for kk = 1:3
            ax = nexttile(tl,kk);
            zmap = nan(size(pmap)); % [rows x shanks]
            for ch = 1:nch
                zmap(chrow(ch),chshank(ch)) = gainvals{con}{kk}(ch);
            end
            h = imagesc(ax, zmap);
            set(h,'AlphaData',isnan(zmap) == false) %no site / removed = blank
            colormap(ax,divmap)
            clim(ax,[-dmax dmax])
            colorbar(ax)
            hold(ax,'on')
            for ch = 1:nch
                if respchans(ch) == true
                    text(ax, chshank(ch), chrow(ch), num2str(chanlist(ch)), 'HorizontalAlignment','center','FontSize',7,'FontWeight','bold')
                else
                    text(ax, chshank(ch), chrow(ch), num2str(chanlist(ch)), 'HorizontalAlignment','center','FontSize',7,'Color',[.4 .4 .4])
                end
            end
            if kk == 3
                sigch = tc_gain_q(:,con) < fdr_q;
                plot(ax, chshank(sigch), chrow(sigch), 'ks', 'MarkerSize',16, 'LineWidth',1.5) %q < fdr_q
            end
            for i = 1:numel(rrow)
                text(ax, rshank(i), rrow(i), 'x', 'HorizontalAlignment','center','Color',[.5 .5 .5])
            end
            hold(ax,'off')
            set(ax,'XTick',1:nshank,'YTick',1:nrow)
            xlabel(ax,'shank')
            ylabel(ax,'row (1 = top)')
            title(ax, gtitles{kk}, 'Interpreter','none')
        end
        title(tl,{[animal ' ' ProbeInfo.Areas{prb} '-Probe: does aligning make the trials more alike than it does on noise? ' condis{con}], ...
            [num2str(sum(resp{con})) ' response trials, window ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms; bold = responsive channel, x = removed; condition ' num2str(con) ' of ' num2str(ncond)]},'Interpreter','none')
    end
    if ncond > 1
        uicontrol(gainfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
            'String','condition:', 'HorizontalAlignment','right')
        uicontrol(gainfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
            'Min',1, 'Max',ncond, 'Value',1, 'SliderStep',[1 1]/(ncond-1), ...
            'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
            num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
    end
    drawnow %show it now; the script doesn't wait, it carries on
    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'tc_gain']);
        pans = findobj(gainfig,'Tag','condpanel');
        for k = 1:ncond
            set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
            pk = findobj(pans,'flat','UserData',k);
            tlk = findobj(pk,'Type','tiledlayout');
            exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' condis{k} '.png'],'Resolution',png_res)
        end
        set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
        set(findobj(gainfig,'Style','slider'),'Value',1);
    end
    %}

    %picture 1: the gain in trial-to-trial r from aligning, aligned conditions only, one box per condition (pairs side
    %by side), one dot per responsive channel, lines join the same channel in the pair's two conditions.
    %Bracket = paired signrank of group 2 vs group 1 over the channels (q = Benjamini-Hochberg over the pairs)
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = axes();
    hold(ax,'on')
    nrc = sum(respchans); %channels (dots) per box
    gall = tc_gainch(respchans,:);
    gtop = max(gall,[],'all','omitnan');
    gbot = min(gall,[],'all','omitnan');
    yr = gtop - gbot; %data range, for spacing the brackets
    if isempty(yr) || isnan(yr) || yr == 0
        yr = 0.01;
        gtop = 0;
        gbot = 0;
    end
    xt = [];
    xtl = {};
    for k = 1:npair
        xpos = (k-1)*3 + [1 2]; %group 1 condition, group 2 condition (a gap between pairs)
        gv = {tc_gainch(respchans,pairs(k,1)), tc_gainch(respchans,pairs(k,2))}; %{2} [responsive channels x 1]
        plot(ax, xpos', [gv{1} gv{2}]', 'Color',[.75 .75 .75], 'LineWidth',0.5) %the same channel in both
        for j = 1:2
            con = pairs(k,j);
            boxchart(ax, xpos(j)*ones(nrc,1), gv{j}, 'BoxFaceColor',condcol(con,:), 'MarkerStyle','none', 'BoxWidth',0.5)
            scatter(ax, xpos(j) + (rand(nrc,1)-0.5)*0.25, gv{j}, 18, condcol(con,:), 'filled') %each channel
            xt(end+1) = xpos(j);
            xtl{end+1} = [condis{con} ' aligned'];
        end
        %the test, as a bracket above the pair
        if isnan(tc_gainpair_p(k)) == false
            yb = max([gv{1}; gv{2}],[],'omitnan') + 0.08*yr; %bracket height
            plot(ax, [xpos(1) xpos(1) xpos(2) xpos(2)], [yb-0.03*yr yb yb yb-0.03*yr], 'k-', 'LineWidth',1)
            if tc_gainpair_q(k) < 0.001
                stars = ' ***';
            elseif tc_gainpair_q(k) < 0.01
                stars = ' **';
            elseif tc_gainpair_q(k) < fdr_q
                stars = ' *';
            else
                stars = '';
            end
            text(ax, mean(xpos), yb + 0.04*yr, sprintf('p = %.3g, q = %.3g%s', tc_gainpair_p(k), tc_gainpair_q(k), stars), ...
                'HorizontalAlignment','center', 'VerticalAlignment','bottom')
        end
    end
    yline(ax,0,':') %no gain
    hold(ax,'off')
    ylim(ax,[gbot - 0.05*yr, gtop + 0.3*yr]) %room for the brackets
    set(ax,'XTick',xt,'XTickLabel',xtl,'TickLabelInterpreter','none')
    ylabel(ax,'gain in trial-to-trial r from aligning (aligned - unaligned mean pairwise r, per channel)')
    title(ax,{[animal ' ' ProbeInfo.Areas{prb} '-Probe: how much aligning the trials raises their consistency (typical-trial template, ' num2str(align_use_ms(1)) '-' num2str(align_use_ms(2)) ' ms)'], ...
        ['one dot per responsive channel (n = ' num2str(nrc) '), lines join the same channel; ' grpnames{2} ' vs ' grpnames{1} ' per pair: paired signrank over the channels, q = Benjamini-Hochberg over the ' num2str(npair) ' pairs']},'Interpreter','none')
    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'tc_gain']);
        exportgraphics(gcf,[save_directory '\Figures_StochRes\' animal '\' fname '.png'],'Resolution',png_res)
    end

    %commented out: one dot per TRIAL (each trial's r averaged over the channels); Roy wants one dot per channel, below
    %{
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
        trval = cell(1,4); %each trial's value: unaligned, aligned, null unaligned, null aligned
        trval{1} = mean(tc_trial{vun}{con}(:,respchans),2,'omitnan'); % [response trials x 1]
        trval{2} = mean(tc_trial{vstr}{con}(:,respchans),2,'omitnan');
        trval{3} = mean(tc_trial_null{vun}{con}(:,respchans),2,'omitnan');
        trval{4} = mean(tc_trial_null{vstr}{con}(:,respchans),2,'omitnan');
        xpos = (j-1)*5 + floor((j-1)/2)*2 + [1 2 3 4]; %unaligned, aligned, null unaligned, null aligned
        ntr = numel(trval{1});
        if ntr == 0
            continue
        end
        %lines joining each trial's unaligned and aligned value (real and null)
        plot(ax, xpos(1:2)', [trval{1} trval{2}]', 'Color',[.75 .75 .75], 'LineWidth',0.5)
        plot(ax, xpos(3:4)', [trval{3} trval{4}]', 'Color',[.85 .85 .85], 'LineWidth',0.5)
        for kk = 1:4
            if kk == 1
                col = condcol(con,:)*0.4 + [1 1 1]*0.6; %lighter for unaligned
            elseif kk == 2
                col = condcol(con,:);
            else
                col = [.6 .6 .6]; %null in grey
            end
            boxchart(ax, xpos(kk)*ones(ntr,1), trval{kk}, 'BoxFaceColor',col, 'MarkerStyle','none', 'BoxWidth',0.5)
            scatter(ax, xpos(kk) + (rand(ntr,1)-0.5)*0.25, trval{kk}, 12, col, 'filled') %each trial
        end
        xt = [xt xpos];
        xtl = [xtl {[condis{con} ' unaligned'], [condis{con} ' aligned'], 'null unal.', 'null al.'}];
    end
    hold(ax,'off')
    set(ax,'XTick',xt,'XTickLabel',xtl,'TickLabelInterpreter','none')
    ylabel(ax,'trial''s mean r with the other trials of its condition (mean over responsive channels)')
    title(ax,{[animal ' ' ProbeInfo.Areas{prb} '-Probe: per-trial consistency, unaligned vs aligned (typical-trial template), real window vs pre-stim null (grey)'], ...
        ['conditions side by side in their pairs (' strjoin(pairnames,', ') '); lines join the same trial before and after aligning']},'Interpreter','none')
    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'tc_trials']);
        exportgraphics(gcf,[save_directory '\Figures_StochRes\' animal '\' fname '.png'],'Resolution',png_res)
    end
    %}

    %picture 3: each responsive CHANNEL's consistency (its mean pairwise trial r within the condition), unaligned vs
    %aligned, real window in the condition's colour, pre-stim null in grey; one dot per channel, lines join the same channel
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    ax = axes();
    hold(ax,'on')
    xt = [];
    xtl = {};
    nrc = sum(respchans); %channels (dots) per box
    porder = reshape(pairs',1,[]); %conditions in pair order: pair 1 A, pair 1 B, pair 2 A, ...
    for j = 1:ncond
        con = porder(j);
        chval = {tc_cond{vun}(respchans,con), tc_cond{vstr}(respchans,con), ...
            tc_cond_null{vun}(respchans,con), tc_cond_null{vstr}(respchans,con)}; %{unaligned, aligned, null unaligned, null aligned} [responsive channels x 1]
        xpos = (j-1)*5 + floor((j-1)/2)*2 + [1 2 3 4]; %unaligned, aligned, null unaligned, null aligned
        if all(isnan(chval{1})) %too few response trials
            continue
        end
        %lines joining each channel's unaligned and aligned value (real and null)
        plot(ax, xpos(1:2)', [chval{1} chval{2}]', 'Color',[.75 .75 .75], 'LineWidth',0.5)
        plot(ax, xpos(3:4)', [chval{3} chval{4}]', 'Color',[.85 .85 .85], 'LineWidth',0.5)
        for kk = 1:4
            if kk == 1
                col = condcol(con,:)*0.4 + [1 1 1]*0.6; %lighter for unaligned
            elseif kk == 2
                col = condcol(con,:);
            else
                col = [.6 .6 .6]; %null in grey
            end
            boxchart(ax, xpos(kk)*ones(nrc,1), chval{kk}, 'BoxFaceColor',col, 'MarkerStyle','none', 'BoxWidth',0.5)
            scatter(ax, xpos(kk) + (rand(nrc,1)-0.5)*0.25, chval{kk}, 18, col, 'filled') %each channel
        end
        xt = [xt xpos];
        xtl = [xtl {[condis{con} ' unaligned'], [condis{con} ' aligned (' num2str(sum(tc_gain_q(respchans,con) < fdr_q)) ' sig)'], 'null unal.', 'null al.'}];
    end
    hold(ax,'off')
    set(ax,'XTick',xt,'XTickLabel',xtl,'TickLabelInterpreter','none')
    ylabel(ax,'channel''s mean pairwise trial r within the condition')
    title(ax,{[animal ' ' ProbeInfo.Areas{prb} '-Probe: trial-to-trial consistency per channel (n = ' num2str(nrc) ' responsive channels per box), unaligned vs aligned (typical-trial template), pre-stim null in grey'], ...
        ['conditions side by side in their pairs (' strjoin(pairnames,', ') '); lines join the same channel; "sig" = channels where aligning beats the null, q < ' num2str(fdr_q)]},'Interpreter','none')
    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'tc_channels']);
        exportgraphics(gcf,[save_directory '\Figures_StochRes\' animal '\' fname '.png'],'Resolution',png_res)
    end
end




%Now beginning the V1 - S1 section



%% Data Analysis: V1 spike-count spectra and S1 LFP spectra (spec_ms)

%one spectrum per channel for the whole spec_ms window (the sustained response), coup_groups conditions only:
%  V1: multitaper spectrum of each response trial's RAW 1 ms spike counts (not the smoothed rate: the 25 ms Gaussian
%      cuts ~half the amplitude by 40 Hz and would bias the low gamma band), divided by what random (Poisson) firing
%      at that trial's own rate gives -> 1 = no rhythm, > 1 = rhythm; mean over the trials
%  S1: multitaper LFP power of the same trials and window, mean over the trials, dB
%Unshifted trials: the same spec_ms window after the real stimulus on every trial. The alignment is only for the trial
%consistency: a trial's own spectrum hardly changes with a shift of a few tens of ms, and the V1 - S1 coupling is
%measured within each trial at the same times on both probes, so it doesn't need lining the trials up either
%(the first draft moved both probes by the trial's lag; dropped, it only slid the window by up to max_lag_ms).
%The tapered spectra of the band frequencies are kept for the coherence in the next section
if do_S1 == 1
    disp('S1: V1 spike-count spectra and S1 LFP spectra')
    tic

    %conditions: coup_groups, with enough response trials for the trial shuffles later
    coupcons = [];
    for con = find(ismember(condgrp, coup_groups))
        if sum(resp{con}) < 3
            run_notes(end+1) = string(condis{con}) + ": fewer than 3 response trials, no V1 - S1 spectra or coupling";
            continue
        end
        coupcons(end+1) = con;
    end
    if isempty(coupcons)
        error('No coup_groups condition has 3 or more response trials, nothing to compare with S1')
    end

    %window: spec_ms, trimmed to the loaded MUA and LFP if it runs past them
    % spec_use_ms = [max([spec_ms(1), load_ms(1)+max_lag_ms, lfp_load_ms(1)+max_lag_ms]), min([spec_ms(2), load_ms(2)-max_lag_ms, lfp_load_ms(2)-max_lag_ms])]; %commented out: the trials aren't shifted here any more, no max_lag_ms room needed
    spec_use_ms = [max([spec_ms(1), load_ms(1), lfp_load_ms(1)]), min([spec_ms(2), load_ms(2), lfp_load_ms(2)])]; %ms
    if isequal(spec_use_ms, spec_ms) == false
        run_notes(end+1) = "Spectrum / coupling window trimmed to fit the loaded data: " + spec_use_ms(1) + " to " + spec_use_ms(2) + " ms (spec_ms was " + spec_ms(1) + " to " + spec_ms(2) + ")";
    end
    sidx = find(t_bin >= spec_use_ms(1) & t_bin < spec_use_ms(2)); %MUA bins of the window
    lidx = find(t_lfp >= spec_use_ms(1) & t_lfp < spec_use_ms(2)); %LFP samples of the window (the same ms, 1 kHz)
    nwin = numel(sidx); %samples per window
    fs1 = 1000; %Hz, the counts and the LFP

    %multitaper tapers, made once (Slepian sequences, unit energy), and the frequencies kept
    ntap = 2*psd_nw - 1;
    tapers = dpss(nwin, psd_nw, ntap); % [window samples x tapers]
    nfft = 2^nextpow2(nwin);
    f_all = (0:nfft/2)' * fs1/nfft; %Hz
    fkeep = f_all >= psd_hz(1) & f_all <= psd_hz(2);
    spec_f = f_all(fkeep); % [freqs x 1] Hz
    nf = numel(spec_f);
    linemask = false(nf,1); %line noise bins, set to NaN if line_out = 1
    if line_out == 1
        for b = 1:size(line_hz,1)
            linemask(spec_f >= line_hz(b,1) & spec_f <= line_hz(b,2)) = true;
        end
    end
    nband = size(bands,1);
    fband = cell(1,nband); %{bands} places in f_all of each band's frequencies
    for b = 1:nband
        fband{b} = find(f_all >= bands(b,1) & f_all <= bands(b,2));
    end
    %(one-sided power spectral density of the windowed stretch, averaged over the tapers: same scaling as pmtm,
    %so random firing at mean count m per ms sits at 2*m/fs)

    rch = find(respchans); %responsive V1 channels (places in chanlist): the V1 side of the grid and the coupling
    nrch = numel(rch);

    v1psd = cell(1,ncond); %{conditions} [freqs x V1 channels] 1 = random firing
    s1psd = cell(1,ncond); %{conditions} [freqs x S1 channels] dB
    Yv = cell(1,ncond); %{conditions}{bands} [band freqs x tapers x trials x responsive V1 channels] tapered spectra, single complex
    Ys = cell(1,ncond); %{conditions}{bands} [band freqs x tapers x trials x S1 channels]
    for con = coupcons
        rtrs = find(resp{con} == true); %response trials (lfp_s1 holds these, same order)
        nresp = numel(rtrs);
        % lb = round(lags{1}{con} / bin_ms); %commented out (and the shifted reads below): no alignment for the spectra / coupling

        %the window of every trial, both probes
        % cnt_al = zeros(nresp,nwin,nch,'single');
        % lfp_al = zeros(nresp,nwin,nch_s1,'single');
        % for tr = 1:nresp
        %     cnt_al(tr,:,:) = spkcount{con}(rtrs(tr), sidx + lb(tr), :);
        %     lfp_al(tr,:,:) = lfp_s1{con}(tr, lidx + lb(tr), :);
        % end
        cnt_al = spkcount{con}(rtrs, sidx, :); % [response trials x window samples x V1 channels] spikes per ms
        lfp_al = lfp_s1{con}(:, lidx, :); % [response trials x window samples x S1 channels] uV

        Yv{con} = cell(1,nband);
        Ys{con} = cell(1,nband);
        for b = 1:nband
            Yv{con}{b} = zeros(numel(fband{b}), ntap, nresp, nrch, 'single');
            Ys{con}{b} = zeros(numel(fband{b}), ntap, nresp, nch_s1, 'single');
        end

        %V1, every kept channel (the responsive ones also keep their band spectra)
        v1psd{con} = nan(nf,nch);
        for ch = 1:nch
            X = double(reshape(cnt_al(:,:,ch), nresp, nwin))'; % [window samples x trials] spikes per ms
            m = mean(X,1); % [1 x trials] mean count per ms of each trial
            P = zeros(nfft/2+1, nresp);
            iv = find(rch == ch); %place among the responsive channels (empty = not responsive)
            for j = 1:ntap
                Y = fft((X - m) .* tapers(:,j), nfft);
                P = P + abs(Y(1:nfft/2+1,:)).^2;
                if isempty(iv) == false
                    for b = 1:nband
                        Yv{con}{b}(:,j,:,iv) = single(reshape(Y(fband{b},:), numel(fband{b}), 1, nresp));
                    end
                end
            end
            P = P / ntap / fs1;
            P(2:end-1,:) = 2*P(2:end-1,:); %one-sided
            lvl = 2*m/fs1; %random-firing level of each trial
            lvl(lvl == 0) = NaN; %a trial without spikes has no spectrum
            R = P(fkeep,:) ./ lvl; % [freqs x trials] 1 = random firing
            v1psd{con}(:,ch) = mean(R,2,'omitnan');
        end

        %S1, every kept channel
        s1psd{con} = nan(nf,nch_s1);
        for c1 = 1:nch_s1
            X = double(reshape(lfp_al(:,:,c1), nresp, nwin))'; % [window samples x trials] uV
            P = zeros(nfft/2+1, nresp);
            for j = 1:ntap
                Y = fft((X - mean(X,1)) .* tapers(:,j), nfft);
                P = P + abs(Y(1:nfft/2+1,:)).^2;
                for b = 1:nband
                    Ys{con}{b}(:,j,:,c1) = single(reshape(Y(fband{b},:), numel(fband{b}), 1, nresp));
                end
            end
            P = P / ntap / fs1;
            P(2:end-1,:) = 2*P(2:end-1,:);
            s1psd{con}(:,c1) = 10*log10(mean(P(fkeep,:),2)); %dB of the trial-mean power
        end
        v1psd{con}(linemask,:) = NaN; %line noise out (if asked)
        s1psd{con}(linemask,:) = NaN;
    end

    %S1: each channel's 1/f background (knee model), so the grid can show the rhythm ABOVE it (0 dB = none): the S1 version
    %of V1's "above random firing". The same fit as NNx_LFPCorrelation_LvsLW.m: log10 power = b - log10(k + f^exponent),
    %first on every point, then again on only the points at or under that first fit (so the rhythm bumps don't drag it
    %up). Fit on each channel's trial-mean spectrum over psd_hz; line noise (line_hz) is always left out of the fit
    notline = true(nf,1); %frequencies the background is fit on
    for b = 1:size(line_hz,1)
        notline(spec_f >= line_hz(b,1) & spec_f <= line_hz(b,2)) = false;
    end
    fitf = spec_f(notline); % [fit freqs x 1] Hz
    apmodel = @(q,f) q(1) - log10(q(2)^2 + f.^q(3)); %knee 1/f model, log10 power; q = [b sqrt(k) exponent], k = q(2)^2 >= 0
    fitopts = optimset('Display','off','MaxFunEvals',3000,'MaxIter',3000);
    apfit1 = @(lp,q0) fminsearch(@(q) sum((lp - apmodel(q,fitf)).^2), q0, fitopts); %fit every point
    aprefit = @(lp,q1) fminsearch(@(q) sum((lp - apmodel(q,fitf)).^2 .* (lp - apmodel(q1,fitf) <= 0)), q1, fitopts); %refit on the points on or under the first fit q1
    apstart = @(c,lp) [lp(end) + log10(8^(-c(1)) + fitf(end)^(-c(1))), sqrt(8^(-c(1))), -c(1)]; %start from a straight line (c = polyfit of log-log): knee at 8 Hz, top end lined up
    s1apq = cell(1,ncond); %{conditions} [S1 channels x 3] b, sqrt(k), exponent
    s1bg = cell(1,ncond); %{conditions} [freqs x S1 channels] dB, the 1/f background
    s1osc = cell(1,ncond); %{conditions} [freqs x S1 channels] dB above the background (0 = no rhythm)
    for con = coupcons
        s1apq{con} = nan(nch_s1,3);
        s1bg{con} = nan(nf,nch_s1);
        s1osc{con} = nan(nf,nch_s1);
        for c1 = 1:nch_s1
            lp = s1psd{con}(notline,c1) / 10; %log10 power (s1psd is dB)
            if all(isfinite(lp)) == false %flat / missing channel stays NaN
                continue
            end
            q = aprefit(lp, apfit1(lp, apstart(polyfit(log10(fitf),lp,1),lp)));
            s1apq{con}(c1,:) = q;
            s1bg{con}(:,c1) = 10*apmodel(q,spec_f); %dB
            s1osc{con}(:,c1) = s1psd{con}(:,c1) - s1bg{con}(:,c1);
        end
    end
    clear cnt_al lfp_al X P Y R
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: V1 spike spectra x S1 LFP spectra grid

%rows = responsive V1 channels (shank by shank, top to bottom), columns = EVERY kept S1 channel (shank by shank, top to
%bottom). Each cell: that V1 channel's spike-count spectrum (black, 1 = random firing, dotted line) and that S1 channel's
%LFP spectrum above its 1/f background (red, 0 dB = dotted line), psd_hz left to right, the bands shaded. One panel per
%condition, slider.
%Drawn as ONE axes per condition, every cell placed in its own box of a grid, so the whole grid is a handful of
%graphics objects (one line for all V1 curves, one for all S1 curves, one patch per band) and draws in seconds
%(was one small axes per pair, thousands of objects, minutes to draw and save)
if do_S1 == 1 && plot_psdgrid == 1
    disp('S1: Making the V1 x S1 spectra grid')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    [~,vord] = sortrows([chshank(rch) chrow(rch)]); %responsive V1 channels shank by shank, top to bottom
    vrows = rch(vord); %places in chanlist, in row order
    [~,scols] = sortrows([chshank_s1 chrow_s1]); %every S1 channel shank by shank, top to bottom (places in chanlist_s1)
    ncol = numel(scols);

    %one y range for all cells and conditions (so they compare): V1 over the responsive channels, S1 over all its channels
    v1all = [];
    s1all = [];
    for con = coupcons
        v1all = [v1all; v1psd{con}(:,rch)];
        s1all = [s1all; s1osc{con}];
    end
    v1lim = [min(v1all,[],'all','omitnan') max(v1all,[],'all','omitnan')];
    s1lim = [min(s1all,[],'all','omitnan') max(s1all,[],'all','omitnan')];
    if any(isnan(v1lim)) || v1lim(1) == v1lim(2)
        v1lim = [0 2];
    end
    if any(isnan(s1lim)) || s1lim(1) == s1lim(2)
        s1lim = [0 1];
    end
    bandcol = [0.3 0.6 1; 1 0.75 0.2]; %shading of band 1, band 2 (more bands reuse these)

    %where everything goes: cell (row r, column c) is the box x = c-1 .. c, y = nrch-r .. nrch-r+1 (row 1 on top),
    %with a small margin inside so neighbouring cells don't touch
    mrg = 0.06; %margin, fraction of a cell
    fx = (spec_f - psd_hz(1)) / diff(psd_hz) * (1 - 2*mrg) + mrg; % [freqs x 1] x of each frequency inside a cell (0..1)
    v1y = @(v) (v - v1lim(1)) / diff(v1lim) * (1 - 2*mrg) + mrg; %V1 value -> height inside a cell (0..1)
    s1y = @(v) (v - s1lim(1)) / diff(s1lim) * (1 - 2*mrg) + mrg; %S1 dB -> height inside a cell
    nf1 = numel(spec_f) + 1; %points per cell in the long lines (+1 for the NaN gap)
    vbreak = find(diff(chshank(vrows)) ~= 0); %rows after which a new V1 shank starts
    sbreak = find(diff(chshank_s1(scols)) ~= 0); %columns after which a new S1 shank starts

    %the reference lines (V1 = 1, S1 = 0 dB) and band patches are the same for every condition: built once
    xref = nan(nf1, nrch*ncol); % [points x cells] cell's x for its lines
    yref1 = nan(2, nrch*ncol); %flat V1 = 1 line per cell (2 points)
    yref0 = nan(2, nrch*ncol); %flat S1 = 0 dB line per cell
    xflat = nan(3, nrch*ncol); %x of the flat lines (2 points + NaN)
    k = 0;
    for r = 1:nrch
        for c = 1:ncol
            k = k + 1;
            xref(1:end-1,k) = (c-1) + fx;
            xflat(1:2,k) = (c-1) + [mrg; 1-mrg];
            yref1(:,k) = (nrch-r) + v1y(1);
            yref0(:,k) = (nrch-r) + s1y(0);
        end
    end
    yref1(3,:) = NaN;
    yref0(3,:) = NaN;
    bandx = cell(1,nband); %{bands} [4 x columns] patch corners: one strip per column, all rows high
    for b = 1:nband
        x0 = (bands(b,1) - psd_hz(1)) / diff(psd_hz) * (1 - 2*mrg) + mrg;
        x1 = (bands(b,2) - psd_hz(1)) / diff(psd_hz) * (1 - 2*mrg) + mrg;
        bandx{b} = (0:ncol-1) + [x0; x1; x1; x0];
    end

    gridfig = figure();
    set(gridfig, 'Position', get(0, 'Screensize'));
    for ic = 1:numel(coupcons)
        con = coupcons(ic);
        pan = uipanel(gridfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
            'Tag','condpanel', 'UserData',ic, 'Visible', ic == 1);
        ax = axes(pan, 'Position',[0.06 0.04 0.92 0.84]);
        hold(ax,'on')

        %band shading, one patch per band
        for b = 1:nband
            patch(ax, bandx{b}, repmat([0; 0; nrch; nrch], 1, ncol), bandcol(mod(b-1,2)+1,:), 'FaceAlpha',0.2, 'EdgeColor','none')
        end

        %every cell's curves, as one long line per probe (cells separated by NaN)
        yv = nan(nf1, nrch*ncol); % [points x cells] V1 curve heights
        ys = nan(nf1, nrch*ncol); %S1 curve heights
        k = 0;
        for r = 1:nrch
            for c = 1:ncol
                k = k + 1;
                yv(1:end-1,k) = (nrch-r) + v1y(v1psd{con}(:,vrows(r)));
                ys(1:end-1,k) = (nrch-r) + s1y(s1osc{con}(:,scols(c)));
            end
        end
        plot(ax, xflat(:), yref1(:), 'k:', 'LineWidth',0.3) %V1 = 1 (random firing)
        plot(ax, xflat(:), yref0(:), 'r:', 'LineWidth',0.3) %S1 = 0 dB (the 1/f background)
        plot(ax, xref(:), yv(:), 'k-', 'LineWidth',0.5)
        plot(ax, xref(:), ys(:), 'r-', 'LineWidth',0.5)

        %cell borders (light), shank borders (dark)
        xline(ax, 0:ncol, 'Color',[.85 .85 .85])
        yline(ax, 0:nrch, 'Color',[.85 .85 .85])
        xline(ax, sbreak, 'Color',[.3 .3 .3], 'LineWidth',1)
        yline(ax, nrch - vbreak, 'Color',[.3 .3 .3], 'LineWidth',1)
        hold(ax,'off')

        xlim(ax,[0 ncol])
        ylim(ax,[0 nrch])
        set(ax, 'XAxisLocation','top', 'XTick',(1:ncol) - 0.5, 'XTickLabel',"CH" + string(chanlist_s1(scols)) + " sh" + string(chshank_s1(scols)), ...
            'YTick',(1:nrch) - 0.5, 'YTickLabel',flipud("V1 CH" + string(chanlist(vrows))), 'TickLength',[0 0], 'FontSize',5, 'Box','on')
        xtickangle(ax,90)
        title(ax,{[animal ' ' condis{con} ': V1 spike-count spectrum (black, dotted = 1 = random firing) vs S1 LFP spectrum above its 1/f background (red, dotted = 0 dB; knee-model fit), all S1 channels (columns, shank by shank)'], ...
            ['each cell: ' num2str(psd_hz(1)) '-' num2str(psd_hz(2)) ' Hz left to right; height V1 ' sprintf('%.2f',v1lim(1)) '-' sprintf('%.2f',v1lim(2)) ' x random firing, S1 ' sprintf('%.1f',s1lim(1)) ' to ' sprintf('%.1f',s1lim(2)) ' dB above 1/f; ' ...
            num2str(sum(resp{con})) ' response trials, ' num2str(spec_use_ms(1)) '-' num2str(spec_use_ms(2)) ' ms; shaded: ' strjoin(bandnames,', ') ' (' char(strjoin((string(bands(:,1)) + "-" + string(bands(:,2)))',', ')) ' Hz); condition ' num2str(ic) ' of ' num2str(numel(coupcons))]},'Interpreter','none','FontSize',8)
    end

    %the slider (finds the panels by their Tag, so it still works in the saved figure)
    if numel(coupcons) > 1
        uicontrol(gridfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
            'String','condition:', 'HorizontalAlignment','right')
        uicontrol(gridfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
            'Min',1, 'Max',numel(coupcons), 'Value',1, 'SliderStep',[1 1]/(numel(coupcons)-1), ...
            'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
            num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
    end
    drawnow %show it now; the script doesn't wait, it carries on

    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'psdgrid']);
        %one PNG per panel (= condition), at twice png_res: the cells are small
        pans = findobj(gridfig,'Tag','condpanel');
        for k = 1:numel(coupcons)
            set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
            pk = findobj(pans,'flat','UserData',k);
            exportgraphics(findobj(pk,'Type','axes'),[save_directory '\Figures_StochRes\' animal '\' fname '-' condis{coupcons(k)} '.png'],'Resolution',2*png_res)
        end
        set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
        set(findobj(gridfig,'Style','slider'),'Value',1);
    end
end




%% Data Analysis: V1 spikes x S1 LFP coupling in the bands, every channel pair

%The SR idea: if the whisker input works as "noise" in a band, V1 spiking and the S1 LFP should go together in that
%band. For every responsive V1 channel x every kept S1 channel, per band (alpha, low gamma), per condition, on the
%response trials in spec_use_ms, unshifted (the same times on both probes; the alignment is only for the trial
%consistency, coupling is measured within each trial):
%  1) spike-field coherence: |cross-spectrum| / sqrt(V1 power x S1 power), averaged over trials and tapers, then over
%     the band's frequencies. Phase coupling: does V1 spike at a consistent S1 phase? (0 = no, 1 = perfect)
%  2) envelope correlation: band-pass both (zero-phase Butterworth), Hilbert amplitude every env_step_ms, per trial
%     r of the two envelopes, Fisher-z mean over trials. Power coupling: does V1 rhythm grow when S1 power grows?
%     (a trial without V1 spikes in the window is left out for that V1 channel)
%  3) (do_mi = 1) Gaussian-copula mutual information (Ince et al. 2017) between the two band signals, each as
%     [real imaginary] (phase and amplitude together, catches non-linear coupling too), samples pooled over time and
%     trials, bias-corrected, bits
%Null: n_shuf trial shuffles (V1 of trial i with S1 of another trial, the same shuffles for every pair and measure).
%Each pair gets a z against its shuffles. Correction over the many pairs: a pair counts when its z beats the
%(1 - pair_p) percentile of the BIGGEST shuffle z over all pairs (per condition, band and measure). Coherence and MI
%are one-sided (more than chance); envelope r is two-sided (|z|).
%Speed: every V1 trial x S1 trial pairing is worked out once per pair (a table), so a shuffle just picks other
%entries of that table
if do_S1 == 1
    disp('S1: V1 spikes x S1 LFP coupling, every channel pair (this takes a while)')
    tic

    nband = size(bands,1);
    esm = sidx(1:env_step_ms/bin_ms:end); %MUA bins of the envelope samples (unshifted)
    esl = lidx(1:env_step_ms:end); %LFP samples of the envelope samples (the same ms)
    nE = numel(esm);
    filt_sos = cell(1,nband); %{bands} band-pass filters (4th order Butterworth, run forwards and backwards)
    filt_g = cell(1,nband);
    for b = 1:nband
        [z,p,k] = butter(2, bands(b,:)/(fs1/2), 'bandpass'); %(order 2 band-pass = 4th order)
        [filt_sos{b},filt_g{b}] = zp2sos(z,p,k);
    end

    %worker pool (Parallel Computing Toolbox only; otherwise nw = 0 = one core, a normal loop in order)
    nw = 0;
    if n_workers > 0 && license('test','Distrib_Computing_Toolbox') && isempty(ver('parallel')) == false
        pool = gcp('nocreate'); %use the pool if one is open
        if isempty(pool)
            pool = parpool; %default size (first time in a session takes a while)
        end
        nw = min(n_workers, pool.NumWorkers);
    end

    coh = cell(1,ncond); %{conditions} [responsive V1 ch x S1 ch x bands]
    envr = cell(1,ncond);
    mi = cell(1,ncond);
    coh_z = cell(1,ncond);
    envr_z = cell(1,ncond);
    mi_z = cell(1,ncond);
    coh_sig = cell(1,ncond);
    envr_sig = cell(1,ncond);
    mi_sig = cell(1,ncond);
    coup_thr = nan(3,nband,ncond); %max-statistic thresholds: coh; envr; mi

    for con = coupcons
        rtrs = find(resp{con} == true); %response trials (lfp_s1 holds these, same order)
        nresp = numel(rtrs);
        % lb = round(lags{1}{con} / bin_ms); %commented out (and the lb(tr) shifts below): coupling uses the unshifted trials
        disp(['S1: ' condis{con} ' (' num2str(nresp) ' response trials, ' num2str(nrch) ' x ' num2str(nch_s1) ' channel pairs)'])

        %the shuffles: V1 of trial i with S1 of trial permset(s,i), no trial with itself; the same for every pair
        rng('default') %same shuffles every run (per condition)
        permset = zeros(n_shuf, nresp);
        for s = 1:n_shuf
            pm = randperm(nresp);
            while any(pm == 1:nresp)
                pm = randperm(nresp);
            end
            permset(s,:) = pm;
        end
        diagidx = (1:nresp)' + ((1:nresp)'-1)*nresp; %table entries of each trial with itself
        shufidx = (1:nresp) + (permset-1)*nresp; % [shuffles x trials] table entries of trial i with its shuffled partner

        %band signals (analytic: amplitude and phase) at the envelope samples, every trial, both probes. The whole loaded
        %trial is filtered, then the aligned window is taken, so the filter edges are far from it
        Av = cell(1,nband); %{bands} [envelope samples x trials x responsive V1 channels] complex, V1 spike counts
        As = cell(1,nband); %{bands} [envelope samples x trials x S1 channels] complex, S1 LFP
        for b = 1:nband
            Av{b} = zeros(nE,nresp,nrch);
            As{b} = zeros(nE,nresp,nch_s1);
        end
        vok = false(nresp,nrch); %V1 trial has spikes in its window (envelope r and MI only use these)
        for iv = 1:nrch
            X = double(reshape(spkcount{con}(rtrs,:,rch(iv)), nresp, nbins))'; % [bins x trials] spikes per ms, whole loaded trial
            % for tr = 1:nresp
            %     vok(tr,iv) = sum(X(sidx + lb(tr), tr)) > 0;
            % end
            vok(:,iv) = sum(X(sidx,:),1)' > 0;
            for b = 1:nband
                H = hilbert(filtfilt(filt_sos{b}, filt_g{b}, X - mean(X,1))); % [bins x trials] analytic band signal
                % for tr = 1:nresp
                %     Av{b}(:,tr,iv) = H(esm + lb(tr), tr);
                % end
                Av{b}(:,:,iv) = H(esm,:);
            end
        end
        for c1 = 1:nch_s1
            X = double(reshape(lfp_s1{con}(:,:,c1), nresp, numel(t_lfp)))'; % [samples x trials] uV, whole loaded trial
            for b = 1:nband
                H = hilbert(filtfilt(filt_sos{b}, filt_g{b}, X - mean(X,1)));
                % for tr = 1:nresp
                %     As{b}(:,tr,c1) = H(esl + lb(tr), tr);
                % end
                As{b}(:,:,c1) = H(esl,:);
            end
        end

        %V1 sides of the tables, made once: envelopes centred and scaled to length 1 per trial (NaN = no spikes)
        Vn = cell(1,nband); %{bands} [envelope samples x trials x responsive V1 channels]
        for b = 1:nband
            E = abs(Av{b});
            E = E - mean(E,1);
            E = E ./ sqrt(sum(E.^2,1));
            for iv = 1:nrch
                E(:, vok(:,iv) == false, iv) = NaN;
            end
            Vn{b} = E;
        end
        %(do_mi) each variable to a standard normal by its rank (the copula), over all samples of all trials together
        Vc = cell(1,nband); %{bands} [envelope samples x trials x 2 (real, imaginary) x responsive V1 channels]
        Sc = cell(1,nband); %{bands} [envelope samples x trials x 2 x S1 channels]
        if do_mi == 1
            for b = 1:nband
                Vc{b} = zeros(nE,nresp,2,nrch);
                Sc{b} = zeros(nE,nresp,2,nch_s1);
                for iv = 1:nrch
                    parts = {real(Av{b}(:,:,iv)), imag(Av{b}(:,:,iv))};
                    for d = 1:2
                        x = parts{d}(:);
                        Vc{b}(:,:,d,iv) = reshape(norminv(tiedrank(x) / (numel(x) + 1)), nE, nresp);
                    end
                end
                for c1 = 1:nch_s1
                    parts = {real(As{b}(:,:,c1)), imag(As{b}(:,:,c1))};
                    for d = 1:2
                        x = parts{d}(:);
                        Sc{b}(:,:,d,c1) = reshape(norminv(tiedrank(x) / (numel(x) + 1)), nE, nresp);
                    end
                end
            end
        end
        Yvc = Yv{con}; %{bands} tapered V1 band spectra (from the spectrum section)
        Ysc = Ys{con};

        %per S1 channel (one per worker): real value + every shuffle's value of every pair and band
        r_coh = nan(nrch,nch_s1,nband); %real values
        r_env = nan(nrch,nch_s1,nband);
        r_mi = nan(nrch,nch_s1,nband);
        s_coh = nan(nrch,nch_s1,nband,n_shuf,'single'); %shuffle values (single: ~n_shuf x pairs x bands x 4 bytes)
        s_env = nan(nrch,nch_s1,nband,n_shuf,'single');
        s_mi = nan(nrch,nch_s1,nband,n_shuf,'single');
        cname = condis{con}; %for the progress lines
        % for c1 = 1:nch_s1 %commented out: on nw workers now (nw = 0 runs it as this loop)
        parfor (c1 = 1:nch_s1, nw) %S1 channels don't depend on each other (each worker gets the V1 tables, a few MB)
            rc = nan(nrch,nband); %this S1 channel's results first, into the big arrays in one go at the end (parfor needs it that way)
            re = nan(nrch,nband);
            rm = nan(nrch,nband);
            sc = nan(nrch,nband,n_shuf,'single');
            se = nan(nrch,nband,n_shuf,'single');
            sm = nan(nrch,nband,n_shuf,'single');
            for b = 1:nband
                B = Ysc{b}(:,:,:,c1); % [band freqs x tapers x trials] S1 tapered spectra
                nfb = size(B,1);
                Syy = mean(abs(B).^2,[2 3]); % [band freqs x 1] S1 power
                Bp = permute(conj(B),[2 3 1]); % [tapers x trials x band freqs]
                Se = abs(As{b}(:,:,c1)); % [envelope samples x trials] S1 envelope, centred and scaled per trial
                Se = Se - mean(Se,1);
                Se = Se ./ sqrt(sum(Se.^2,1));
                for iv = 1:nrch
                    %1) coherence: table of every V1 trial x S1 trial cross-spectrum (summed over tapers)
                    A = Yvc{b}(:,:,:,iv); % [band freqs x tapers x trials]
                    Sxx = mean(abs(A).^2,[2 3]); % [band freqs x 1] V1 power
                    T = reshape(pagemtimes(permute(A,[3 2 1]), Bp), nresp*nresp, nfb) / ntap; % [(trials*trials) x band freqs]
                    den = sqrt(Sxx .* Syy)'; % [1 x band freqs]
                    rc(iv,b) = mean(abs(mean(T(diagidx,:),1)) ./ den);
                    Tsh = reshape(T(shufidx(:),:), n_shuf, nresp, nfb); %every shuffle's pairings
                    sc(iv,b,:) = single(mean(abs(reshape(mean(Tsh,2), n_shuf, nfb)) ./ den, 2));

                    %2) envelope correlation: table of every V1 trial x S1 trial r
                    R = Vn{b}(:,:,iv)' * Se; % [trials x trials] (NaN rows = V1 trials without spikes)
                    R = max(min(R, 0.9999), -0.9999); %keep atanh finite
                    Z = atanh(R);
                    re(iv,b) = tanh(mean(Z(diagidx),'omitnan'));
                    se(iv,b,:) = single(mean(Z(shufidx),2,'omitnan')); %Fisher z means (the real one is compared as z too, below)

                    %3) mutual information, Gaussian copula, V1 [re im] vs S1 [re im], V1 trials with spikes only
                    if do_mi == 1
                        vi = find(vok(:,iv)); %V1 trials used
                        nv = numel(vi);
                        N = nE*nv; %samples
                        Va = reshape(Vc{b}(:,vi,:,iv), N, 2); % [samples x 2] V1, fixed for every shuffle
                        Cxx = cov(Va);
                        mx = mean(Va,1)';
                        %per-trial pieces of the S1 side and of the cross products, so any pairing is just sums
                        Sall = Sc{b}(:,:,:,c1); % [envelope samples x trials x 2]
                        Gs = zeros(2,2,nresp); %each S1 trial's sum of products
                        ms = zeros(2,nresp); %each S1 trial's sums
                        for tr = 1:nresp
                            st = reshape(Sall(:,tr,:), nE, 2);
                            Gs(:,:,tr) = st'*st;
                            ms(:,tr) = sum(st,1)';
                        end
                        Txy = zeros(nresp,nresp,2,2); %V1 trial i x S1 trial j sum of products, per dimension pair
                        for d1 = 1:2
                            for d2 = 1:2
                                Txy(:,:,d1,d2) = Vc{b}(:,:,d1,iv)' * Sall(:,:,d2);
                            end
                        end
                        %bias correction for Gaussian entropies from N samples (Ince et al.): the same for real and shuffles
                        psiterms = psi((N - (1:4))/2) / 2;
                        bcorr = psiterms(3) + psiterms(4) - psiterms(1) - psiterms(2);
                        partners = [vi'; permset(:,vi)]; % [1 + shuffles x V1 trials] S1 partner of each V1 trial (row 1 = real)
                        mival = zeros(n_shuf+1,1);
                        for s = 1:n_shuf+1
                            pj = partners(s,:);
                            my = sum(ms(:,pj),2) / N;
                            Cyy = (sum(Gs(:,:,pj),3) - N*(my*my')) / (N-1);
                            Sxy = zeros(2,2);
                            for d1 = 1:2
                                for d2 = 1:2
                                    tt = Txy(:,:,d1,d2);
                                    Sxy(d1,d2) = sum(tt(vi' + (pj-1)*nresp));
                                end
                            end
                            Cxy = (Sxy - N*(mx*my')) / (N-1);
                            Cj = [Cxx Cxy; Cxy' Cyy];
                            mival(s) = (0.5*log(det(Cxx)*det(Cyy)/det(Cj)) + bcorr) / log(2); %bits
                        end
                        rm(iv,b) = mival(1);
                        sm(iv,b,:) = single(mival(2:end));
                    end
                end
            end
            r_coh(:,c1,:) = reshape(rc, nrch, 1, nband);
            r_env(:,c1,:) = reshape(re, nrch, 1, nband);
            r_mi(:,c1,:) = reshape(rm, nrch, 1, nband);
            s_coh(:,c1,:,:) = reshape(sc, nrch, 1, nband, n_shuf);
            s_env(:,c1,:,:) = reshape(se, nrch, 1, nband, n_shuf);
            s_mi(:,c1,:,:) = reshape(sm, nrch, 1, nband, n_shuf);
            disp(['S1: ' cname ', coupling shuffles done for S1 CH' num2str(chanlist_s1(c1))])
        end

        %z against the shuffles, and the max-over-pairs threshold (per band and measure)
        %(envelope r is tested on the Fisher-z mean, the same scale as its shuffles)
        realv = {r_coh, atanh(max(min(r_env,0.9999),-0.9999)), r_mi};
        shufv = {s_coh, s_env, s_mi};
        zout = cell(1,3);
        sigout = cell(1,3);
        for meas = 1:3
            if meas == 3 && do_mi == 0
                continue
            end
            sh = double(shufv{meas}); % [V1 ch x S1 ch x bands x shuffles]
            mu = mean(sh,4,'omitnan');
            sd = std(sh,0,4,'omitnan');
            zp = (realv{meas} - mu) ./ sd; % [V1 ch x S1 ch x bands]
            zs = (sh - mu) ./ sd; %every shuffle as z
            sg = false(size(zp));
            for b = 1:nband
                zsb = reshape(zs(:,:,b,:), [], n_shuf); % [pairs x shuffles]
                if meas == 2 %two-sided
                    maxz = max(abs(zsb),[],1,'omitnan');
                    thr = prctile(maxz, 100*(1 - pair_p));
                    sg(:,:,b) = abs(zp(:,:,b)) > thr;
                else %one-sided
                    maxz = max(zsb,[],1,'omitnan');
                    thr = prctile(maxz, 100*(1 - pair_p));
                    sg(:,:,b) = zp(:,:,b) > thr;
                end
                coup_thr(meas,b,con) = thr;
            end
            zout{meas} = zp;
            sigout{meas} = sg;
        end
        coh{con} = r_coh;
        envr{con} = r_env;
        coh_z{con} = zout{1};
        envr_z{con} = zout{2};
        coh_sig{con} = sigout{1};
        envr_sig{con} = sigout{2};
        if do_mi == 1
            mi{con} = r_mi;
            mi_z{con} = zout{3};
            mi_sig{con} = sigout{3};
        end
    end
    clear Av As Vn Vc Sc Yvc Ysc X H E s_coh s_env s_mi sh zs realv shufv
    if mr_clean == 1
        clear Yv Ys lfp_s1 %the tapered spectra and the S1 LFP aren't needed past here
    end
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end



%% Data Visualization: V1 x S1 coupling heatmaps

%one panel per condition (slider): a heatmap per measure (rows) x band (columns), rows = responsive V1 channels,
%columns = S1 channels (both shank by shank, top to bottom, thin lines between shanks). Colour = the measure itself
%(coherence, envelope r, MI in bits; one colour range per measure x band shared by all conditions), black squares =
%significant against the trial shuffles after the max-over-pairs correction.
%NOTE: the raw values include the chance level (few trials and sparse V1 channels push coherence / MI up, the shared
%stimulus response pushes envelope r up), so a strong colour alone isn't coupling: trust the squares (z is saved: *_z)
%(was coloured by z vs the shuffles; Roy wanted the measure itself)
%With both groups in coup_groups, extra panels: LW - L of the measure itself (not z), per position-matched pair
if do_S1 == 1 && plot_coup == 1
    disp('S1: Making the coupling heatmaps')

    if plotsave == 1
        %create figure folder
        if isfolder([save_directory '\Figures_StochRes\' animal]) == false
            mkdir([save_directory '\Figures_StochRes\' animal])
        end
    end

    nc = 128; %blue - white - red colour map (0 = white)
    divmap = [linspace(0,1,nc)' linspace(0.3,1,nc)' ones(nc,1); ones(nc,1) linspace(1,0.2,nc)' linspace(1,0.2,nc)'];
    [~,vord] = sortrows([chshank(rch) chrow(rch)]); %responsive V1 channels shank by shank, top to bottom (rows of the results are in rch order)
    [~,sord] = sortrows([chshank_s1 chrow_s1]); %S1 channels shank by shank, top to bottom
    vbreak = find(diff(chshank(rch(vord))) ~= 0) + 0.5; %between V1 shanks
    sbreak = find(diff(chshank_s1(sord)) ~= 0) + 0.5; %between S1 shanks
    measnames = {'spike-field coherence','envelope correlation','mutual information'};
    nmeas = 2 + do_mi;

    %panels: each condition, then (both groups) each pair with both conditions analysed
    panlist = zeros(0,2); % [panels x 2] condition (or pair), 1 = condition / 2 = difference
    for con = coupcons
        panlist(end+1,:) = [con 1];
    end
    for k = 1:npair
        if all(ismember(pairs(k,:), coupcons))
            panlist(end+1,:) = [k 2];
        end
    end
    npan = size(panlist,1);
    ptags = cell(1,npan); %file name tag of each panel

    %one colour range per measure x band over all the conditions (so the slider compares like with like):
    %coherence and MI from 0 up, envelope r symmetric around 0
    vv = {coh, envr, mi}; %{measures}{conditions} [V1 ch x S1 ch x bands]
    rawlim = ones(nmeas,nband); % [measures x bands] top of the colour range
    for meas = 1:nmeas
        for b = 1:nband
            allv = [];
            for con = coupcons
                allv = [allv; reshape(vv{meas}{con}(:,:,b),[],1)];
            end
            if meas == 2
                rawlim(meas,b) = max(abs(allv),[],'omitnan');
            else
                rawlim(meas,b) = max(allv,[],'omitnan');
            end
            if isempty(rawlim(meas,b)) || isnan(rawlim(meas,b)) || rawlim(meas,b) <= 0
                rawlim(meas,b) = 1;
            end
        end
    end
    cbunits = {'coherence','envelope r (Fisher-z mean over trials)','MI (bits)'};

    coupfig = figure();
    set(coupfig, 'Position', get(0, 'Screensize'));
    for ip = 1:npan
        pan = uipanel(coupfig, 'Units','normalized', 'Position',[0 0.05 1 0.95], 'BorderType','none', ...
            'Tag','condpanel', 'UserData',ip, 'Visible', ip == 1);
        tl = tiledlayout(pan, nmeas, nband, 'TileSpacing','compact','Padding','compact');
        for meas = 1:nmeas
            for b = 1:nband
                ax = nexttile(tl,(meas-1)*nband + b);
                if panlist(ip,2) == 1 %one condition: the measure itself
                    con = panlist(ip,1);
                    % zz = {coh_z{con}, envr_z{con}, mi_z{con}}; %commented out (and the z lines below): colour = the measure now
                    ss = {coh_sig{con}, envr_sig{con}, mi_sig{con}};
                    % M = zz{meas}(vord,sord,b);
                    M = vv{meas}{con}(vord,sord,b); % [V1 rows x S1 columns]
                    S = ss{meas}(vord,sord,b);
                    % lim = max([max(abs(M),[],'all','omitnan') coup_thr(meas,b,con)]);
                    % cbtxt = 'z vs trial shuffles';
                    if meas == 2 %envelope r: either sign
                        cmapuse = divmap;
                        climuse = [-rawlim(meas,b) rawlim(meas,b)];
                    else %coherence, MI: 0 and up
                        cmapuse = parula;
                        climuse = [0 rawlim(meas,b)];
                    end
                    cbtxt = cbunits{meas};
                    ttxt = [measnames{meas} ', ' bandnames{b} ' (' num2str(bands(b,1)) '-' num2str(bands(b,2)) ' Hz): ' num2str(sum(S,'all')) ' of ' num2str(numel(S)) ' pairs significant (thr z ' sprintf('%.2f',coup_thr(meas,b,con)) ')'];
                else %pair difference, group 2 - group 1, of the measure itself
                    ca = pairs(panlist(ip,1),1);
                    cb = pairs(panlist(ip,1),2);
                    M = vv{meas}{cb}(vord,sord,b) - vv{meas}{ca}(vord,sord,b);
                    S = false(size(M));
                    lim = max(abs(M),[],'all','omitnan');
                    if isempty(lim) || isnan(lim) || lim == 0
                        lim = 1;
                    end
                    cmapuse = divmap;
                    climuse = [-lim lim];
                    cbtxt = [condis{cb} ' - ' condis{ca}];
                    ttxt = [measnames{meas} ', ' bandnames{b} ' (' num2str(bands(b,1)) '-' num2str(bands(b,2)) ' Hz): ' condis{cb} ' - ' condis{ca}];
                end
                % if isempty(lim) || isnan(lim) || lim == 0 %commented out: each branch sets its own colour range now
                %     lim = 1;
                % end
                h = imagesc(ax, M);
                set(h,'AlphaData',isnan(M) == false)
                % colormap(ax,divmap)
                % clim(ax,[-lim lim])
                colormap(ax,cmapuse)
                clim(ax,climuse)
                cbar = colorbar(ax);
                cbar.Label.String = cbtxt;
                hold(ax,'on')
                [sr,sc] = find(S);
                plot(ax, sc, sr, 's', 'MarkerEdgeColor','k', 'MarkerSize',5, 'LineWidth',1) %significant pairs
                for i = 1:numel(vbreak)
                    yline(ax, vbreak(i), 'Color',[.3 .3 .3])
                end
                for i = 1:numel(sbreak)
                    xline(ax, sbreak(i), 'Color',[.3 .3 .3])
                end
                hold(ax,'off')
                set(ax,'YTick',1:nrch,'YTickLabel',string(chanlist(rch(vord))),'XTick',1:nch_s1,'XTickLabel',string(chanlist_s1(sord)),'FontSize',5)
                ylabel(ax,'V1 channel (responsive)','FontSize',8)
                xlabel(ax,'S1 channel','FontSize',8)
                title(ax, ttxt, 'Interpreter','none','FontSize',8)
            end
        end
        if panlist(ip,2) == 1
            con = panlist(ip,1);
            title(tl,{[animal ' ' condis{con} ': V1 spikes x S1 LFP coupling, every channel pair'], ...
                [num2str(sum(resp{con})) ' response trials, ' num2str(spec_use_ms(1)) '-' num2str(spec_use_ms(2)) ' ms, ' num2str(n_shuf) ' trial shuffles, max-over-pairs p < ' num2str(pair_p) ' (black squares); panel ' num2str(ip) ' of ' num2str(npan)]},'Interpreter','none')
            ptags{ip} = condis{con};
        else
            title(tl,{[animal ' ' pairnames{panlist(ip,1)} ': V1 spikes x S1 LFP coupling, group 2 - group 1 (the measures, not z; no test)'], ...
                ['panel ' num2str(ip) ' of ' num2str(npan)]},'Interpreter','none')
            ptags{ip} = strrep(pairnames{panlist(ip,1)},' ','_');
        end
    end
    if npan > 1
        uicontrol(coupfig, 'Style','text', 'Units','normalized', 'Position',[0.18 0.008 0.11 0.03], ...
            'String','panel:', 'HorizontalAlignment','right')
        uicontrol(coupfig, 'Style','slider', 'Units','normalized', 'Position',[0.3 0.01 0.4 0.03], ...
            'Min',1, 'Max',npan, 'Value',1, 'SliderStep',[1 1]/(npan-1), ...
            'Callback', @(src,evt) set(findobj(src.Parent,'Tag','condpanel'), {'Visible'}, ...
            num2cell(cell2mat(get(findobj(src.Parent,'Tag','condpanel'),'UserData')) == round(src.Value))));
    end
    drawnow %show it now; the script doesn't wait, it carries on

    if plotsave == 1
        fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' 'coupling']);
        pans = findobj(coupfig,'Tag','condpanel');
        for k = 1:npan
            set(pans, {'Visible'}, num2cell([pans.UserData]' == k));
            pk = findobj(pans,'flat','UserData',k);
            tlk = findobj(pk,'Type','tiledlayout');
            exportgraphics(tlk,[save_directory '\Figures_StochRes\' animal '\' fname '-' ptags{k} '.png'],'Resolution',png_res)
        end
        set(pans, {'Visible'}, num2cell([pans.UserData]' == 1));
        set(findobj(coupfig,'Style','slider'),'Value',1);
    end
end




%Now beginning the report and saving


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
    if resp_method == 2
        report(end+1) = "Response trials (method 2: channels 1-" + resp_shank_nch + " of one shank all stay >= " + resp_z + " robust SDs above their all-trial baseline " + base_use_ms(1) + " to " + base_use_ms(2) + " ms for >= " + resp_min_ms + " ms, within " + resp_together_ms + " ms of each other, " + peak_win_ms + " ms windows, in " + resp_ms(1) + "-" + resp_ms(2) + " ms):";
    else
        report(end+1) = "Response trials (template matching, " + resp_ms(1) + "-" + resp_ms(2) + " ms; score > the " + resp_null_prct + "th percentile of the condition's fake no-stim scores):";
    end
    for con = 1:ncond
        report(end+1) = "  " + condis{con} + ": " + sum(resp{con}) + " of " + ntrials(con) + " (" + sprintf('%.0f',100*respfrac(con)) + "%), threshold " + sprintf('%.2f',respthr(con)) + ...
            ", " + sum(respchans_cond(:,con)) + " responsive channels, " + tmpl_rounds(con) + " template round(s)";
        %the spike-count rule
        if exist('respcount','var')
            dcnt = respcount{con}(resp_detect{con}); %counts of the detected trials
            if isempty(dcnt)
                cntxt = "no detected trials";
            else
                cntxt = "detected trials' count median " + sprintf('%.0f',median(dcnt)) + ", 10th percentile " + sprintf('%.0f',prctile(dcnt,10));
            end
            report(end+1) = "    detected " + sum(resp_detect{con}) + ", dropped " + resp_lowcount(con) + " for < " + resp_mincount + " spikes (" + resp_min_spk + "/channel x " + sum(respchans) + " responsive channels) in " + resp_count_ms(1) + "-" + resp_count_ms(2) + " ms; " + cntxt;
        end
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
if exist('lagsd','var')
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
if exist('tc_cond','var')
    report(end+1) = "Trial-to-trial consistency, mean pairwise r over the responsive channels (pre-stim null in brackets), window " + align_use_ms(1) + "-" + align_use_ms(2) + " ms:";
    for con = 1:ncond
        txt = "  " + condis{con} + ":";
        for v = 1:nver
            txt = txt + "  " + al_names{v} + " " + sprintf('%.3f',mean(tc_cond{v}(respchans,con),'omitnan')) + " (" + sprintf('%.3f',mean(tc_cond_null{v}(respchans,con),'omitnan')) + ")";
        end
        txt = txt + ";  aligning beats the null (q < " + fdr_q + ") on " + sum(tc_gain_q(:,con) < fdr_q) + " of " + nch + " channels (" + sum(tc_gain_q(respchans,con) < fdr_q) + " of " + sum(respchans) + " responsive)";
        report(end+1) = txt;
    end
    report(end+1) = "  condition vs condition (label shuffles, Benjamini-Hochberg over channels x comparisons), channels with q < " + fdr_q + " (median r difference over the responsive channels):";
    for i = 1:ncmp
        txt = "    " + cmpnames{i} + ":";
        for v = 1:nver
            dd = median(tc_cond{v}(respchans,cmppairs(i,2)) - tc_cond{v}(respchans,cmppairs(i,1)),'omitnan');
            txt = txt + "  " + al_names{v} + " " + sum(tc_cmp_q{v}(:,i) < fdr_q) + " (" + sprintf('%+.3f',dd) + ")";
        end
        report(end+1) = txt;
    end
    report(end+1) = "  gain in r from aligning (aligned - unaligned, per responsive channel), " + grpnames{2} + " vs " + grpnames{1} + " per pair (paired signrank over the channels, Benjamini-Hochberg over the pairs):";
    for k = 1:npair
        report(end+1) = "    " + pairnames{k} + ": median gain " + sprintf('%.4f',median(tc_gainch(respchans,pairs(k,1)),'omitnan')) + " vs " + sprintf('%.4f',median(tc_gainch(respchans,pairs(k,2)),'omitnan')) + ", median difference " + sprintf('%+.4f',tc_gainpair_d(k)) + ", p = " + sprintf('%.3g',tc_gainpair_p(k)) + ", q = " + sprintf('%.3g',tc_gainpair_q(k)) + ", " + sum(respchans) + " channels";
    end
end

%V1 x S1 coupling
if do_S1 == 1 && exist('coh_sig','var')
    report(end+1) = "V1 spikes x S1 LFP coupling (" + nrch + " responsive V1 x " + nch_s1 + " S1 channels, " + spec_use_ms(1) + "-" + spec_use_ms(2) + " ms, unshifted trials, " + n_shuf + " trial shuffles, max-over-pairs p < " + pair_p + "):";
    measnames = {'coherence','envelope r','MI (bits)'};
    for con = coupcons
        report(end+1) = "  " + condis{con} + " (" + sum(resp{con}) + " response trials):";
        kneef = (s1apq{con}(:,2).^2) .^ (1 ./ s1apq{con}(:,3)); %Hz, knee frequency of each S1 channel
        report(end+1) = "    S1 1/f background over the S1 channels: median exponent " + sprintf('%.2f',median(s1apq{con}(:,3),'omitnan')) + ", median knee " + sprintf('%.1f',median(kneef,'omitnan')) + " Hz";
        vals = {coh{con}, envr{con}, mi{con}};
        sigs = {coh_sig{con}, envr_sig{con}, mi_sig{con}};
        zs = {coh_z{con}, envr_z{con}, mi_z{con}};
        for meas = 1:2 + do_mi
            for b = 1:size(bands,1)
                sg = sigs{meas}(:,:,b);
                vb = vals{meas}(:,:,b);
                txt = "    " + measnames{meas} + ", " + bandnames{b} + ": median " + sprintf('%.3f',median(vb,'all','omitnan')) + ", " + sum(sg,'all') + " of " + numel(sg) + " pairs significant";
                if any(sg,'all')
                    zb = abs(zs{meas}(:,:,b));
                    zb(sg == false) = NaN;
                    [~,ip] = max(zb(:)); %strongest significant pair
                    [iv,c1] = ind2sub(size(zb), ip);
                    txt = txt + "; strongest V1 CH" + chanlist(rch(iv)) + " x S1 CH" + chanlist_s1(c1) + " (" + sprintf('%.3f',vb(ip)) + ", z " + sprintf('%.1f',zs{meas}(iv,c1,b)) + ")";
                end
                report(end+1) = txt;
            end
        end
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
    'stim_onset_30k','stim_onset_lfp','fs_mua','win_ms','base_ms','pad_ms','bin_ms','smooth_ms', ...
    'resp_ms','pre_resp_ms','pre_end_ms','fake_ms','tmpl_trim','tmpl_iter','chan_z_min','resp_null_prct','resp_method','base_use_ms','peak_win_ms','resp_z','resp_shank_nch','resp_min_ms','resp_together_ms','resp_count_ms','resp_min_spk','keep_qc','qc_plotsave', ...
    'align_ms','max_lag_ms','n_perm','fdr_q', ...
    'do_S1','s1_region','lfp_load_ms','coup_groups','spec_ms','psd_nw','psd_hz','line_out','line_hz', ...
    'bands','bandnames','env_step_ms','do_mi','n_shuf','pair_p','load_ms'};
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
    'align_methods','align_use_ms','lags','lags_null','lagpk','lagsd','lagsd_null','t_plot','rate_mean_al','rate_sem_al', ...
    'al_names','tc_trial','tc_trial_null','tc_cond','tc_cond_null','cmppairs','cmpnames','tc_cmp_p','tc_cmp_q','tc_p', ...
    'tc_gain','tc_gain_p','tc_gain_q','tc_alvsun_p','tc_gainch','tc_gainpair_p','tc_gainpair_q','tc_gainpair_d', ...
    'prb_s1','pmap_s1','chanlist_s1','chrow_s1','chshank_s1','t_lfp','coupcons','spec_use_ms','spec_f','v1psd','s1psd','s1apq','s1bg','s1osc','rch', ...
    'coh','envr','mi','coh_z','envr_z','mi_z','coh_sig','envr_sig','mi_sig','coup_thr'};
%Quality Control details: how the response trials and responsive channels were found (keep_qc)
if keep_qc == 1
    savelist = [savelist, {'chan_medz','chan_medz_cond','respchans_cond','respscore','respnull','respthr','resp_tmpl','tmpl_rounds', ...
        'align_ref','align_ref_r','chpeakz','resp_shankch','resp_detect','respcount','resp_lowcount','resp_mincount'}];
end
haveit = false(1,numel(savelist));
for i = 1:numel(savelist)
    haveit(i) = exist(savelist{i},'var') == 1;
end
savelist = savelist(haveit);

fname = sprintf([animal '-StochRes_MUA']);
save([save_directory '\StochRes\' animal '\' fname], savelist{:}, '-v7.3');
disp(['saved: ' save_directory '\StochRes\' animal '\' fname '.mat'])

