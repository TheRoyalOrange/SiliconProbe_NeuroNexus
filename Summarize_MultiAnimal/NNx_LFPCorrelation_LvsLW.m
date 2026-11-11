% NNx_LFPCorrelation_LvsLW.m
%
% Description: follow-up to NNx_P1PeakTime_Variance_LvsLW.m (archived 2026-11-08, out of the repo but
%   still in git history), straight from the LFP files. Asks whether
%   adding whisker (LW) makes the response a clearer, more consistent oscillation than light alone (L),
%   three ways, for the channels of one probe (region, e.g. V1):
%     1) autocorrelation (ACF) of each trial in the P2 window (350-3000 ms post stim): mean ACF trace per
%        channel and condition (probe-layout plot, L vs LW overlaid), plus per-trial numbers: 1/e decay,
%        area, height and lag of the first secondary peak (rhythmicity), and how fast the successive
%        peaks die out (damping rate).
%     2) trial-to-trial correlation (1-3000 ms) within each channel and condition: how alike the trials
%        are (stim-locked consistency). Per channel L vs LW by a label-shuffle permutation test.
%     3) channel-to-channel correlation (1-3000 ms) per trial: L, LW and LW - L matrices per animal,
%        and each channel's mean correlation with the other channels.
%     4) spectrum of each trial in the P2 window (Welch, 1 s Hann segments, bandpass gain divided back
%        out): 1/f background fit with the knee model (log10 power = b - log10(k + f^exponent), FOOOF-style
%        refit on the points at or under the first fit, line noise left out), its exponent, and the
%        oscillation power above it (dB) in one band per channel x pair (biggest bump of L + LW pooled,
%        +/- band_hw), so L and LW are measured in the same band. Added 2026-11-08 because single-trial ACF
%        peaks at pk_minprom = 0.02 show up on pure 1/f noise (every trial "has" one) and follow the
%        background slope rather than a rhythm.
%     5) trial classes (added 2026-11-08), worked out on every good trial BEFORE the analyses above:
%        - rhythmic response: the trial's P2 oscillation power (4) beats the same measure on a same-length
%          pre-stim window (rhythm_null_prct of the pre-stim windows, mean z over the oscillating channels)
%        - P1 response: the trial's P1 (P1_ms) matches its channels' average P1 waveform (leave-one-out
%          template) better than p1_null_prct of fake no-stim windows in the pre-stim baseline do
%        With rhythm_filter = 1 everything in 1) to 4) runs on the rhythmic trials only (the question is the
%        STRUCTURE of the rhythmic response, and trial means / fractions otherwise shift with how often the
%        animal responds at all). The classes are compared L vs LW (stats_class, Venn diagrams): rhythmic out
%        of all trials and out of the P1-response trials only.
%   You list the animals and the L and LW condition names in the setup (paired in list order, like
%   the P1 script), then pick each condition's LFP files (all animals at once, matched by file name).
%   Everything is done per animal, then combined across animals like the P1 script: per "unit"
%   (animal x channel x pair) L vs LW, signrank across units and a mixed model with animal and channel
%   as random effects, optional per-animal channel maps, ranked-shank heatmaps
%   (NNx_plotLWIntensityCompare/NNx_plotLWTTXCompare_multianimal) and a per-trial x channel csv for R.
%   ACF secondary peak vs trial-trial correlation: the ACF sees rhythm whatever its phase, the
%   trial-trial correlation only sees what is time-locked to the stimulus, so together they split
%   "LW is more rhythmic" from "LW is more time-locked".
%
% Inputs:
%   User-set config (Setup section): mice, region, condis_L, condis_LW, tr_conditional_use, windows, filter, spectrum
%     (spec_fit_hz, line_hz, osc_hz, band_hw, all Hz), trial classes (rhythm_filter, rhythm_base_end_ms, P1_ms,
%     p1_base_ms, fake_ms (ms), rhythm_null_prct, p1_null_prct (percentiles), chan_z_min (z)), toggles
%   LFP files picked per condition with uipickfiles, <lfp_dir>\<animal>\<animal>-<stim>-LFP.mat:
%     stim_lfp_stimchunks (double, trials x samples x channels) - peri-stimulus LFP, 1 kHz,
%       presumably uV [inferred: units not stated upstream]
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1)
%     tr_remove (0/1, 1 x nTrials) - mask over tr_keep (updated format from
%       NNx_TrialRemove_Conditional(_callable).m); only tr_keep(tr_remove == 0) are used. Old format
%       (different length) stops the script with an error
%     tr_remove_conditional (table, optional) - Name (string) / TrialIdx (double, tr_keep values; NaN =
%       evaluated, none flagged) named trial sets; trials under any tr_conditional_use name are left out too.
%       A file never evaluated for a name keeps all its trials and gets a note in the report
%   <probeinfo_dir>\<animal>-ProbeInfo.mat, ProbeInfo fields:
%     .Areas (cell of char, 1 x probes) - region per probe; region picks the probe
%     .ProbeMaps (cell, 1 x probes) - [rows x shanks] raw channel IDs, 0 = no site, row 1 = top
%     .Ch_Remove (cell, 1 x probes) - raw channel IDs to leave out
%     .ShankLabels (cell, 1 x probes, optional) - struct of per-shank scores -> ShankRank_<score> columns
%   INFERRED DATA CONTRACTS (as in NNx_SummaryAnalysis_MultiMouse_allChannels.m, not stated in the files):
%     - stim onset is sample 5001 of stim_lfp_stimchunks (BaseAnalysis epochs with pre = 5000)
%     - channel index (dim 3) of stim_lfp_stimchunks = raw channel ID (true when ChanIds = 1:Chans)
%     - Ch_Remove is indexed by probe number
%     - the 5000 pre-stim samples are there and the pre-stim baseline (-4600 to 0 ms) is free of the
%       previous trial's response (pre-stim rhythm window and fake P1 windows)
%
% Outputs (workspace; nothing saved unless save_csv / plotsave / save_heatpng):
%   lfp (cell, animals x conditions) - (single, trials x keep_ms x channels) analyzed trials (rhythmic only
%     when rhythm_filter = 1), uV, bandpassed if lfpfilt = 1, files of the same animal + condition pooled
%   lfp_base (cell, animals x conditions) - (single, trials x nP2 x channels) pre-stim window of every good trial
%   trinfo (cell, animals x conditions) - (double, trials x 2) Condition_FileNum, Trial (tr_keep value), analyzed trials
%   ntrials_all, ntrials (double, animals x conditions) - good trials / analyzed trials
%   oscP2, oscBase (cell, animals x conditions) - (double, trials x channels) dB above the
%     trial's own background in the shared band, pre-stim / P2 (oscP2 subset to the analyzed trials)
%   appq_trbase (cell, animals x conditions) - (trials x channels x 3) pre-stim background fits, every good trial
%   p1score, p1resp (cell, animals x conditions) - (trials x 1) P1 template score (mean z over the responsive
%     channels) and true = P1 response, every good trial; p1null, p1chans (cell, 1 x animals), p1thr (1 x animals)
%   rhythmscore, rhythmic (cell, animals x conditions) - (trials x 1) rhythm score (mean pre-stim z over the
%     oscillating channels) and true = rhythmic response, every good trial; rhythmnull, oscchans (cell, animals x
%     pairs), rhythmthr (animals x pairs)
%   trialclass (table, one row per good trial) - Animal_Name, Condition_Name, Condition_FileNum, Trial, LvsLW,
%     Pair, P1_Score, P1_Resp (0/1), Rhythm_Score, Rhythmic (0/1), TrialClass (Rhythmic+P1 / Rhythmic only /
%     P1 only / Neither)
%   stats_class (table) - per pair x animal: L and LW counts (n, rhythmic, P1, both, neither), % rhythmic out of
%     all trials and out of the P1-response trials, Fisher exact p of each L vs LW
%   stats_class_glme (table) - per pair (>= 2 animals): fitglme(Rhythmic ~ Cond + (1|Animal), binomial) on the
%     P1-response trials, odds ratio LW / L with CI and p
%   csv (save_csv = 1) - also trialclass to <csv_directory>\<name>_trialclass.csv
%   acf (cell, animals x conditions) - (double, trials x lags x channels) P2 autocorrelation, 1 at lag 0
%   acfmean, acfsem (cell, animals x pairs x 2 groups [L LW]) - (double, lags x channels) mean ACF / SEM
%   acfmet (cell, animals x conditions) - (double, trials x channels x 5) per-trial ACF metrics, in
%     the order of metricnames(1:5)
%   acfmet_mean (cell, animals x pairs x 2) - (double, channels x 5) the same metrics on the mean ACF
%   trialcorr (cell, animals x conditions) - (double, trials x channels) the trial's mean r with the other
%     trials of its group (Fisher z mean, back to r)
%   tcmean, tc_p (cell, animals x pairs) - (channels x 2) mean pairwise trial r in L / LW; (channels x 1)
%     permutation p of LW - L (two-sided)
%   chcorr (cell, animals x conditions) - (double, trials x channels) channel's mean r with the other
%     channels on that trial; offdiag (trials x 1) mean r over all channel pairs on that trial
%   chanZsum (cell, animals x conditions) - (double, channels x channels) sum over trials of Fisher z
%     channel-channel correlation (divide by the trial count for the mean)
%   fitf (double, freqs x 1) - Hz the 1/f background is fit on (spec_fit_hz minus line_hz, 1 Hz steps)
%   specP (cell, animals x conditions) - (single, trials x fitf x channels) log10 P2 power, bandpass gain
%     divided out
%   specmeanP (cell, animals x pairs x 2) / specpool (cell, animals x pairs) - (double, channels x fitf)
%     log10 of the mean power of the group / of every trial of the pair
%   appq_pool, appq_mean (channels x 3), appq_tr (trials x channels x 3) - knee-model fits, q = [b sqrt(k)
%     exponent]: log10 power = q(1) - log10(q(2)^2 + f^q(3)); knee frequency = (q(2)^2)^(1/q(3)) Hz
%   specband (cell, animals x pairs) - (double, channels x 1) centre of each channel's shared band, Hz
%   specmet (cell, animals x conditions) - (double, trials x channels x 3) exponent, oscillation power
%     (dB), shared band centre (Hz)
%   specmet_mean (cell, animals x pairs x 2) - (double, channels x 3) the same on the group's mean
%     spectrum, except column 3 = the group's own peak frequency (Hz)
%   specres (cell, animals x conditions) - (single, trials x fitf x channels) dB above the trial's own fit
%   linedb (cell, animals x conditions) - (double, trials x channels) line-noise spike over 3-6 Hz
%     either side, dB (report only)
%   tableres (table, trials x channels rows) - one row per good trial x channel, columns:
%     Animal_Name, Condition_Name, Region, LvsLW, Pair (string); Condition_FileNum, Trial, Channel_ID,
%     Chan_Shank, Chan_Depth (1 = top kept channel of the shank), ShankRank_<score> (double), Rhythmic and
%     P1_Resp (0/1, the trial's class; Rhythmic all 1 when rhythm_filter = 1), and
%     P2_ACFDecay (ms) - first lag where the ACF drops below 1/e; NaN if not within acf_maxlag_ms
%     P2_ACFArea (ms) - sum of the ACF over lags 1..acf_maxlag_ms (1 ms steps)
%     P2_ACFPk2 (unitless) - height of the first ACF peak after lag 0 with prominence >= pk_minprom
%       (the first bump after the first trough); NaN if none. Higher = more rhythmic
%     P2_ACFPk2Lag (ms) - lag of that peak (1000/lag = dominant frequency, Hz)
%     P2_ACFPkDamp (1/s) - how fast the positive ACF peaks shrink: minus the slope of log(peak height)
%       vs lag, x1000. Needs >= 2 positive peaks, else NaN. Smaller = longer-lasting oscillation
%     All_TrialCorr (r) - see trialcorr (window corr_ms, named All_ so the heatmap functions take it)
%     All_ChanCorr (r) - see chcorr
%     P2_SpecExp (unitless) - 1/f exponent of the trial's own background fit (above the knee)
%     P2_SpecOscPow (dB) - mean dB above the trial's background in the channel x pair shared band
%     P2_SpecPkFreq (Hz) - per trial: the shared band centre (same for every trial of a channel x pair);
%       unit level (_L/_LW): each group's own peak on its mean spectrum, no per-trial test (_p = NaN)
%   unitres (table, units x columns) - one row per animal x channel x pair: n trials, fraction of trials
%     with a secondary ACF peak, and per metric <metric>_L, <metric>_LW (channel-level: ACF metrics on
%     the mean ACF, mean pairwise r, mean channel r), <metric>_p (ranksum on per-trial values; the
%     permutation test for All_TrialCorr), <metric>_q (Benjamini-Hochberg within animal x pair)
%   stats_paired (table) - per metric x pair (+ "all"): medians, n units LW higher, two-sided signrank
%     across units, and the mixed model fitlme(value ~ Cond (+ Pair) + (1 + Cond|Animal) + (1|ChanKey),
%     REML, Satterthwaite DF) on log (P2_ACFDecay, P2_ACFPk2Lag; effect = LW/L ratio), Fisher z
%     (correlations; effect = z difference) or raw values (the rest; effect = LW - L)
%   lme_models (cell) - the fitted models, same rows as stats_paired
%   stats_animal (table) - per pair x animal x metric: channels, medians, channels LW higher, signrank p
%     across that animal's channels, channels with q < 0.05
%   stats_chanchan (table) - per pair x animal: mean channel-pair r in L / LW, trials, ranksum p over trials
%   run_notes (string) - problems hit along the way, shown in the report
%   (heatplot > 0) heatfigs, heatdata (cell) - from the heatmap functions, one per alignment label
%   csv (save_csv = 1) - tableres to <csv_directory>\<name>.csv
%   figures (plotsave = 1) - <save_directory>\Figures_LFPCorr\<animal>\ and \Figures_LFPCorr\
%
% Dependencies: NNx_BaseAnalysis*.m (LFP files, ProbeInfo); NNx_TrialRemove_Conditional(_callable).m
%   (tr_remove mask format, tr_remove_conditional sets); ProbeInfo.ShankLabels from the channel/shank labelling scripts (heatmaps);
%   NNx_plotLWIntensityCompare_multianimal.m / NNx_plotLWTTXCompare_multianimal.m (heatplot); uipickfiles.m;
%   Signal Processing Toolbox (butter, zp2sos, filtfilt, pwelch, freqz); Statistics Toolbox (corr, ranksum,
%   signrank, fitlme, fishertest, fitglme); MATLAB R2022a+ (clim, tiledlayout 'tight', islocalmax).
%

%% Setup: Specify Animals, Conditions, and Analysis Info

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Who is the data? (exactly how their LFP file names start: <animal>-<stim>-LFP.mat)
mice = {"20260508-p9", "20260509-p10", "20260821-p10"}; %"20260226-p12", "20260411-p9", "20260423-p12" 

%which probe should I analyze? (must match an entry of ProbeInfo.Areas for every mouse, or a probe number like "1")
region = "V1";

%which conditions are LIGHT (L) and which are LIGHT + WHISKER (LW)? (you'll pick each one's files next)
%same number of each -> paired in list order (1st L with 1st LW, ...) and compared pair by pair.
%otherwise all L vs all LW as one lumped pair. Name them like the summary csv (L_8, LW_8, L_8_TTX...)
%if you want the ranked-shank heatmaps
%condis_L = {'L_4','L_8','L_12', 'L_15'};
%condis_LW = {'LW_4','LW_8','LW_12', 'LW_15'};
condis_L = {'L_8','L_8_TTX'};
condis_LW = {'LW_8','LW_8_TTX'};

%leave out trials flagged under these conditional trial set names? (made with NNx_TrialRemove_Conditional.m,
%e.g. {"whisker_twitch_artifact"}). {} = no, only tr_remove
tr_conditional_use = {};

%where are the LFP and ProbeInfo files?
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%Where should I, Matlab, put the results?
csv_directory = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR';
save_directory = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal'; %figures go in \Figures_LFPCorr

%%%Analysis details%%%
stim_onset = 5001; %sample of stimulus onset (0 ms) in stim_lfp_stimchunks, 1 kHz
P2_ms = [350 3000]; %ms post stimulus, window for the autocorrelation
corr_ms = [1 3000]; %ms post stimulus, window for the trial-trial and channel-channel correlations
acf_maxlag_ms = 1000; %longest autocorrelation lag, ms (long enough to see a few secondary peaks)

%bandpass the LFP first? (slow drift puffs up the autocorrelation and hides the secondary peaks)
lfpfilt = 1;
    % 0 = raw LFP (each window is only demeaned)
    % 1 = zero-phase Butterworth bandpass of lfp_band, on a padded stretch so the window edges are clean
lfp_band = [2 100]; %Hz
filt_pad_ms = 1000; %ms of extra data on each side for the filter to settle

%take out the evoked response first? (autocorrelation, spectrum and channel-channel only, never trial-trial)
remove_erp = 0;
    % 0 = keep it (evoked + induced)
    % 1 = subtract each condition's mean trace from its trials (induced activity only)

pk_minprom = 0.02; %how much a bump in the autocorrelation has to stick out (prominence) to count as a peak
n_perm = 1000; %label shuffles for the per-channel trial-trial correlation test

%%%Spectrum details (1/f background + oscillation power above it, P2 window)%%%
spec_fit_hz = [3 90]; %Hz, range the 1/f background is fit over (keep it inside lfp_band if lfpfilt = 1)
line_hz = [48 52]; %Hz, line noise left out of everything spectral (50 Hz mains). One row per band, e.g. [48 52; 98 102]
osc_hz = [3 45]; %Hz, where to look for each channel's oscillation band
band_hw = 2; %Hz, half width of that band (peak +/- band_hw)

%%%Trial classes (rhythmic response? P1 response?)%%%
%only analyze trials with a rhythmic response? (rhythmic = the P2 oscillation power beats the same measure
%on a same-length pre-stim window, see Data Analysis: trial classes)
rhythm_filter = 1;
    % 0 = no, every good trial (the trial classes are still worked out and reported)
    % 1 = yes, ACF / spectrum / trial-trial / channel-channel statistics only on rhythmic trials
rhythm_base_end_ms = -50; %ms, the pre-stim window ends here and is as long as P2_ms (-2700 to -50 ms)
rhythm_null_prct = 95; %a trial is rhythmic when its score beats this percentile of the pre-stim scores (~1 window per trial, so 99 is too noisy)
P1_ms = [75 350]; %ms post stimulus, P1 response window (same as the summary script)
p1_base_ms = 200; %ms just before the (real or fake) stimulus, the P1 is measured from its mean
fake_ms = -4400:350:-400; %ms, fake stimulus times in the pre-stim baseline for the P1 no-stim reference
p1_null_prct = 99; %a trial has a P1 response when its score beats this percentile of the fake no-stim scores
chan_z_min = 1; %a channel counts as oscillating / responding when its typical (median) trial z is at least this

%%%Plotting details%%%
%trial class pictures? (score vs no-stim histograms, Venn of rhythmic vs P1 trials, class counts)
plot_trialclass = 1;
    % 0 = no
    % 1 = yes
%probe-layout autocorrelation traces? (one figure per animal x pair)
plot_acf = 1;
    % 0 = no
    % 1 = yes
%probe-layout spectra above the 1/f background? (one figure per animal x pair)
plot_spec = 1;
    % 0 = no
    % 1 = yes
%per-animal channel maps of LW - L? (one figure per metric x pair, one probe tile per animal, like the P1 script)
plot_chanmaps = 1;
chanmap_metrics = ["P2_ACFPk2","P2_ACFPkDamp","P2_SpecOscPow","All_TrialCorr","All_ChanCorr"]; %any of the metric names in the header
%channel x channel correlation matrices? (one figure per animal x pair)
plot_chanchan = 1;
%across-animal summary (paired L vs LW per unit, one panel per metric)?
plot_summary = 1;
%Want the ranked-shank summary heatmaps too?
heatplot = 2;
    % 0 = no heatmaps
    % 1 = intensity layout (NNx_plotLWIntensityCompare_multianimal: L/LW 4, 8, 12, 15)
    % 2 = TTX layout (NNx_plotLWTTXCompare_multianimal: L/LW 8 with and without TTX)
heat_metrics = ["P2_ACFPk2","P2_SpecOscPow","All_TrialCorr","All_ChanCorr"]; %one heatmap figure each
save_heatpng = 0; %1 = save the heatmaps as PNG (you'll get a file name box)
%save the per-trial table as a csv for R?
save_csv = 0;
    % 0 = no
    % 1 = yes (you'll get a file name box)
%save the figures?
plotsave = 0;

%%%nothing else to fill in%%%


%% Data PreProcessing: line up the conditions into L/LW pairs

mice = string(mice); %string array, easier to compare
nmice = numel(mice);
condis = [condis_L condis_LW]; %all conditions, L first
ncond = numel(condis);
condi_L = "L"; %labels used in the tables and plots
condi_LW = "LW";

if numel(condis_L) == numel(condis_LW)
    pairnames = string(condis_L) + " vs " + string(condis_LW); %one pair per L/LW match
    npairs = numel(pairnames);
    grpconds = cell(npairs,2); %{pair x group} condition indices, group 1 = L, 2 = LW
    for j = 1:npairs
        grpconds{j,1} = j;
        grpconds{j,2} = numel(condis_L) + j;
    end
else
    pairnames = strjoin(string(condis_L),"+") + " vs " + strjoin(string(condis_LW),"+"); %everything lumped together
    npairs = 1;
    grpconds = {1:numel(condis_L), numel(condis_L) + (1:numel(condis_LW))};
end
condpair = zeros(1,ncond); %which pair each condition is in
condgrp = zeros(1,ncond); %and which group (1 = L, 2 = LW)
for j = 1:npairs
    for g = 1:2
        condpair(grpconds{j,g}) = j;
        condgrp(grpconds{j,g}) = g;
    end
end

run_notes = strings(0,1); %problems for the report at the end


%% Data Extraction: pick the LFP files of each condition (all animals at once)

superCondis_dir = cell(1,ncond);
for con = 1:ncond
    superCondis_dir{con} = uipickfiles('FilterSpec',lfp_dir,'Prompt',['Choose ' condis{con} ' LFP files (all animals)'],'REFilter','-LFP\.mat$');
    if iscell(superCondis_dir{con}) == false
        error(['No files picked for ' condis{con}])
    end
end

%which mouse is each file from? (file names start with animal-)
condiFileMice = cell(1,ncond); %per file, index into mice
condiFileNum = cell(1,ncond); %per file, Nth file of that mouse in that condition
for con = 1:ncond
    for file = 1:numel(superCondis_dir{con})
        [~,fname] = fileparts(superCondis_dir{con}{file});
        m = [];
        for i = 1:nmice
            if startsWith(fname, mice(i) + "-")
                m(end+1) = i; %#ok<SAGROW>
            end
        end
        if numel(m) ~= 1
            error(['Condition ' condis{con} ': file "' fname '" matches ' num2str(numel(m)) ' of the listed mice (needs exactly 1)'])
        end
        condiFileMice{con}(file) = m;
        condiFileNum{con}(file) = sum(condiFileMice{con} == m); %Nth file of this mouse in this condition
    end
end


%% Data Extraction: each animal's probe, channels, and shank labels (ProbeInfo)

pinfos = cell(1,nmice); %each animal's ProbeInfo
mouseProbe = zeros(1,nmice); %probe index used per animal
for m = 1:nmice
    pf = fullfile(probeinfo_dir, mice(m) + "-ProbeInfo.mat");
    if isfile(pf) == false
        error(['No ProbeInfo file for ' char(mice(m)) ': ' char(pf)])
    end
    pinfo = load(pf,'ProbeInfo');
    pinfos{m} = pinfo.ProbeInfo;
    prb = find(strcmp(string(pinfos{m}.Areas), region));
    if isempty(prb) %region can also be a probe number
        prb = str2double(region);
    end
    if numel(prb) ~= 1 || isnan(prb)
        error([char(mice(m)) ': region "' char(region) '" does not pick exactly one probe in ProbeInfo.Areas'])
    end
    mouseProbe(m) = prb;
end
clear pinfo

%which shank scores are there? (ShankLabels fields with one number per shank) -> ShankRank_<score> columns
scorenames = strings(1,0);
for m = 1:nmice
    SL = [];
    if isfield(pinfos{m},'ShankLabels') && numel(pinfos{m}.ShankLabels) >= mouseProbe(m)
        SL = pinfos{m}.ShankLabels{mouseProbe(m)};
    end
    if isstruct(SL)
        nshank = size(pinfos{m}.ProbeMaps{mouseProbe(m)},2);
        fn = fieldnames(SL);
        for f = 1:numel(fn)
            if isnumeric(SL.(fn{f})) && numel(SL.(fn{f})) == nshank
                scorenames(end+1) = string(fn{f}); %#ok<SAGROW>
            end
        end
    end
end
scorenames = unique(scorenames,'stable');
shankscore_use = strings(1,0);
if isempty(scorenames) == false
    [sel,ok] = listdlg('ListString',cellstr(scorenames),'SelectionMode','multiple','Name','Shank scores', ...
        'PromptString','Shank scores to add as ShankRank_ columns (for the heatmaps):','ListSize',[350 200], ...
        'InitialValue',1:numel(scorenames));
    if ok == 1
        shankscore_use = scorenames(sel);
    end
end
rankcols = matlab.lang.makeValidName("ShankRank_" + shankscore_use);

%kept channels per animal (probe minus Ch_Remove), where they sit, and their shank ranks
chanlist = cell(1,nmice); %column of kept raw channel IDs (sorted)
chrow = cell(1,nmice); %ProbeMaps row of each kept channel
chshank = cell(1,nmice); %ProbeMaps column (shank)
chdepth = cell(1,nmice); %position among the shank's kept channels, 1 = top
chorder = cell(1,nmice); %channel order shank by shank, top to bottom (for the channel x channel plots)
chrank = cell(1,nmice); %[channels x shank scores] rank of the channel's shank, 1 = highest score
for m = 1:nmice
    prb = mouseProbe(m);
    pmap = pinfos{m}.ProbeMaps{prb}; %[rows x shanks] raw channel IDs, 0 = no site
    chremove = [];
    if numel(pinfos{m}.Ch_Remove) >= prb
        chremove = pinfos{m}.Ch_Remove{prb};
    end
    chanlist{m} = setdiff(nonzeros(pmap), chremove(:)); %sorted column of raw channel IDs

    %kept channels numbered 1..n from the top of each shank (same as the summary script)
    depthmap = nan(size(pmap));
    for s = 1:size(pmap,2)
        kept = pmap(:,s) > 0 & ismember(pmap(:,s), chremove) == false;
        depthmap(kept,s) = 1:nnz(kept);
    end

    nch = numel(chanlist{m});
    chrow{m} = zeros(nch,1);
    chshank{m} = zeros(nch,1);
    for ch = 1:nch
        [chrow{m}(ch), chshank{m}(ch)] = find(pmap == chanlist{m}(ch)); %where this channel sits on the probe
    end
    chdepth{m} = depthmap(sub2ind(size(pmap), chrow{m}, chshank{m}));
    [~,chorder{m}] = sortrows([chshank{m} chrow{m}]);

    %shank ranks: highest score = 1, NaN score -> NaN rank, ties -> lower shank number first
    chrank{m} = nan(nch,numel(shankscore_use));
    SL = [];
    if isfield(pinfos{m},'ShankLabels') && numel(pinfos{m}.ShankLabels) >= prb
        SL = pinfos{m}.ShankLabels{prb};
    end
    for k = 1:numel(shankscore_use)
        sname = char(shankscore_use(k));
        if isstruct(SL) && isfield(SL, sname) && isnumeric(SL.(sname)) && numel(SL.(sname)) == size(pmap,2)
            score = double(SL.(sname));
            shankrank = nan(1,numel(score));
            valid = find(isnan(score) == false);
            [~,ord] = sort(score(valid),'descend');
            shankrank(valid(ord)) = 1:numel(valid);
            chrank{m}(:,k) = shankrank(chshank{m})';
        else
            run_notes(end+1) = mice(m) + ": no ShankLabels." + shankscore_use(k) + ", its ShankRank column is NaN"; %#ok<SAGROW>
        end
    end
end
clear SL


%% Data Extraction: load the good trials of every file

%only the stretch both windows need is kept (plus padding for the filter while filtering)
keep_ms = min(P2_ms(1),corr_ms(1)) : max(P2_ms(2),corr_ms(2)); %ms post stimulus kept in lfp
P2idx = find(keep_ms >= P2_ms(1) & keep_ms <= P2_ms(2)); %P2 window within keep_ms
corridx = find(keep_ms >= corr_ms(1) & keep_ms <= corr_ms(2)); %correlation window within keep_ms

if lfpfilt == 1
    pad = filt_pad_ms; %samples = ms at 1 kHz
    [z,p,k] = butter(2, lfp_band/(1000/2), 'bandpass'); %2nd order per edge, as the summary script
    [lfp_sos,lfp_g] = zp2sos(z,p,k); %second-order sections, stable at low cutoffs
else
    pad = 0;
end
%loadsamp = stim_onset + ((keep_ms(1)-pad):(keep_ms(end)+pad)); %samples pulled from each file %commented out: the trial classes need the pre-stim part too
%inkeep = pad + (1:numel(keep_ms)); %where keep_ms sits inside loadsamp
loadsamp = 1:(stim_onset + keep_ms(end) + pad); %from the start of the trial (pre-stim baseline, for the trial classes) to the end of the windows + padding
t_load = loadsamp - stim_onset; %ms of each loaded sample
inkeep = find(t_load >= keep_ms(1) & t_load <= keep_ms(end)); %where keep_ms sits inside loadsamp
nP2 = numel(P2idx); %samples (ms) in the P2 window
inbase = find(t_load > rhythm_base_end_ms - nP2 & t_load <= rhythm_base_end_ms); %pre-stim window as long as P2 (rhythm reference)
inp1 = find(t_load >= P1_ms(1) & t_load <= P1_ms(2)); %P1 window
inp1base = find(t_load > -p1_base_ms & t_load <= 0); %just before the stimulus
nfake = numel(fake_ms);
%(the earliest fake P1 windows start ~400 ms into the trial, inside the filter's settling stretch; the 2 Hz
%highpass settles in ~0.1 s, so that's fine)

lfp = cell(nmice,ncond); %{animal x condition} [trials x keep_ms x kept channels], uV, single
lfp_base = cell(nmice,ncond); %{animal x condition} [trials x nP2 x kept channels] pre-stim window, uV, single
p1w = cell(nmice,ncond); %{animal x condition} [trials x P1 samples x kept channels] P1 minus the mean just before the stimulus
p1fake = cell(nmice,ncond); %{animal x condition} [(trials*fake windows) x P1 samples x kept channels] the same at the fake times
trinfo = cell(nmice,ncond); %{animal x condition} [trials x 2] Condition_FileNum, Trial
%filesused = strings(0,3); %commented out: 4th column = trials left out by tr_conditional_use
filesused = strings(0,4); %animal, condition, file, trials left out by tr_conditional_use (for the report)

tic
for con = 1:ncond
    for file = 1:numel(superCondis_dir{con})
        m = condiFileMice{con}(file);
        [~,fname] = fileparts(superCondis_dir{con}{file});
        disp(['Loading: ' condis{con} ', ' fname])
        d = load(superCondis_dir{con}{file},'stim_lfp_stimchunks','tr_keep','tr_remove');
        if numel(d.tr_remove) ~= numel(d.tr_keep)
            error([fname ': tr_remove is not a 0/1 mask over tr_keep (old format). Run the file through NNx_TrialRemove_Conditional first'])
        end
        if loadsamp(1) < 1 || loadsamp(end) > size(d.stim_lfp_stimchunks,2)
            error([fname ': the windows + filter padding run off the end of the trials (' num2str(size(d.stim_lfp_stimchunks,2)) ' samples)'])
        end
        trs = d.tr_keep(d.tr_remove == 0); %good trials
        trs = trs(:);
        %leave out the trials flagged under tr_conditional_use (TrialIdx = tr_keep values, NaN = evaluated, none flagged)
        nexcl = 0;
        if isempty(tr_conditional_use) == false
            if any(strcmp(who('-file',superCondis_dir{con}{file}),'tr_remove_conditional'))
                tc = load(superCondis_dir{con}{file},'tr_remove_conditional');
                tc = tc.tr_remove_conditional;
            else
                tc = table('Size',[0,2],'VariableTypes',{'string','double'},'VariableNames',{'Name','TrialIdx'});
            end
            for k = 1:numel(tr_conditional_use)
                namerows = strcmp(tc.Name, tr_conditional_use{k});
                if any(namerows) == false
                    run_notes(end+1) = string(fname) + ": never evaluated for " + tr_conditional_use{k} + ", none of its trials left out"; %#ok<SAGROW>
                    continue
                end
                flagged = ismember(trs, tc.TrialIdx(namerows));
                nexcl = nexcl + sum(flagged);
                trs = trs(flagged == false);
            end
        end
        dat = d.stim_lfp_stimchunks(trs, loadsamp, chanlist{m}); %[good trials x samples x kept channels]
        clear d

        [ntr,nsamp,nch] = size(dat);
        if lfpfilt == 1
            X = reshape(permute(dat,[2 1 3]), nsamp, []); %[samples x (trials*channels)], filtfilt works down columns
            X = filtfilt(lfp_sos, lfp_g, X); %zero phase
            dat = permute(reshape(X, nsamp, ntr, nch), [2 1 3]); %back to [trials x samples x channels]
            clear X
        end

        lfp{m,con} = cat(1, lfp{m,con}, single(dat(:,inkeep,:))); %pool the files of the same animal + condition
        lfp_base{m,con} = cat(1, lfp_base{m,con}, single(dat(:,inbase,:)));
        p1w{m,con} = cat(1, p1w{m,con}, single(dat(:,inp1,:) - mean(dat(:,inp1base,:),2)));
        nf = zeros(ntr,nfake,numel(inp1),nch,'single'); %[trials x fake windows x P1 samples x channels]
        for f = 1:nfake
            fw = t_load >= fake_ms(f) + P1_ms(1) & t_load <= fake_ms(f) + P1_ms(2); %fake P1 window
            fb = t_load > fake_ms(f) - p1_base_ms & t_load <= fake_ms(f); %and just before it
            nf(:,f,:,:) = reshape(dat(:,fw,:) - mean(dat(:,fb,:),2), ntr, 1, numel(inp1), nch);
        end
        p1fake{m,con} = cat(1, p1fake{m,con}, reshape(nf, [], numel(inp1), nch));
        clear nf
        trinfo{m,con} = cat(1, trinfo{m,con}, [repmat(condiFileNum{con}(file),numel(trs),1) trs]);
        filesused(end+1,:) = [mice(m) string(condis{con}) string(superCondis_dir{con}{file}) string(nexcl)]; %#ok<SAGROW>
        clear dat
    end
end
elapsed_time = toc;
fprintf('loading execution time: %.2f seconds\n', elapsed_time);

%trial counts, and animals missing a condition
ntrials = zeros(nmice,ncond); %[animals x conditions]
for m = 1:nmice
    for con = 1:ncond
        ntrials(m,con) = size(lfp{m,con},1);
        if ntrials(m,con) == 0
            run_notes(end+1) = mice(m) + ": no trials for " + string(condis{con}); %#ok<SAGROW>
        end
    end
end

%induced-only version for the autocorrelation and channel-channel parts (remove_erp = 1)
if remove_erp == 1
    lfp_ind = cell(nmice,ncond);
    for m = 1:nmice
        for con = 1:ncond
            if ntrials(m,con) > 0
                lfp_ind{m,con} = lfp{m,con} - mean(lfp{m,con},1); %each trial minus the condition's mean trace
            end
        end
    end
else
    lfp_ind = lfp; %same data (no copy is made unless one gets changed)
end


%% Data Analysis 1c, part 1: spectrum of every trial (P2 + a same-length pre-stim window), 1/f background, oscillation power

%An oscillation is a bump in the power spectrum sitting on the smooth 1/f background. The background gets
%fit with the knee model (like FOOOF: log10 power = b - log10(k + f^exponent)), first on every point, then
%again on only the points at or under that first fit, so the bumps don't drag it up. Oscillation power =
%how far the spectrum sits above the background, in dB. A change in the background's slope goes into the
%exponent instead of faking a rhythm (which is what fools the single-trial ACF peaks).
%Each channel gets ONE band per pair: the biggest bump of L + LW pooled (picked without knowing which trial
%is which), peak +/- band_hw. L and LW trials are both measured in that same band.
%The same is done on a pre-stim window as long as P2 (ending at rhythm_base_end_ms), so each trial's rhythm can be
%held up against what the channel does before any stimulus (Data Analysis: trial classes). This part uses EVERY
%good trial; the group numbers come after the rhythm filter (Data Analysis 1c, part 2)
disp('Spectrum: fitting the 1/f background of every trial and channel (P2 and pre-stim)')
tic
spec_win = 1000; %samples (= ms) per Welch segment -> 1 Hz steps, 4 half-overlapping segments in the P2 window
[~,spec_f] = pwelch(zeros(nP2,1),hann(spec_win),spec_win/2,spec_win,1000); %frequencies, Hz
if lfpfilt == 1
    filtpow = abs(lfp_g * freqz(lfp_sos,spec_f,1000)).^4; %power the bandpass took (filtfilt runs it twice -> |H|^4), divided back out below
else
    filtpow = ones(size(spec_f));
end
notline = true(size(spec_f)); %leave the line noise out
for i = 1:size(line_hz,1)
    notline(spec_f >= line_hz(i,1) & spec_f <= line_hz(i,2)) = false;
end
fsel = find(spec_f >= spec_fit_hz(1) & spec_f <= spec_fit_hz(2) & notline); %rows of spec_f the background is fit on
fitf = spec_f(fsel); %Hz
nfit = numel(fitf);
oscm = fitf >= osc_hz(1) & fitf <= osc_hz(2); %where to look for each channel's band (within fitf)
fosc = fitf(oscm);
iline = find(abs(spec_f - mean(line_hz(1,:))) < 0.5); %line frequency bin (50 Hz)
nbline = abs(spec_f - spec_f(iline)) >= 3 & abs(spec_f - spec_f(iline)) <= 6; %its neighbours, 3-6 Hz either side

%every trial's spectrum, P2 and pre-stim
specP = cell(nmice,ncond); %{animal x condition} [trials x fitf x channels] log10 P2 power, filter gain divided out
specPbase = cell(nmice,ncond); %{animal x condition} [trials x fitf x channels] the same for the pre-stim window
linedb = cell(nmice,ncond); %{animal x condition} [trials x channels] line-noise spike over its neighbours, dB
for m = 1:nmice
    nch = numel(chanlist{m});
    for con = 1:ncond
        if ntrials(m,con) == 0
            continue
        end
        specP{m,con} = zeros(ntrials(m,con),nfit,nch,'single');
        specPbase{m,con} = zeros(ntrials(m,con),nfit,nch,'single');
        linedb{m,con} = zeros(ntrials(m,con),nch);
        for ch = 1:nch
            x = double(lfp_ind{m,con}(:,P2idx,ch))'; %[P2 samples x trials], pwelch works down columns
            x = x - mean(x,1); %demean each trial's window
            Pw = pwelch(x,hann(spec_win),spec_win/2,spec_win,1000) ./ filtpow; %[freqs x trials]
            specP{m,con}(:,:,ch) = log10(Pw(fsel,:))';
            linedb{m,con}(:,ch) = 10*log10(Pw(iline,:) ./ mean(Pw(nbline,:),1))';
            x = double(lfp_base{m,con}(:,:,ch))'; %same for the pre-stim window
            x = x - mean(x,1);
            Pw = pwelch(x,hann(spec_win),spec_win/2,spec_win,1000) ./ filtpow;
            specPbase{m,con}(:,:,ch) = log10(Pw(fsel,:))';
        end
    end
end

%L + LW pooled spectrum of each channel (every trial of the pair), for the shared band. Mean in power, then log
specpool = cell(nmice,npairs); %{animal x pair} [channels x fitf] log10 power
for m = 1:nmice
    for j = 1:npairs
        if sum(ntrials(m,grpconds{j,1})) == 0 || sum(ntrials(m,grpconds{j,2})) == 0 %pair is missing L or LW
            continue
        end
        S = cat(1, specP{m,[grpconds{j,:}]}); %[trials of the pair x fitf x channels]
        specpool{m,j} = log10(reshape(mean(10.^double(S),1),nfit,[]))';
    end
end
clear S

%the background fit, written once so every fit uses the very same steps. q = [b sqrt(k) exponent], k = q(2)^2 >= 0
apmodel = @(q,f) q(1) - log10(q(2)^2 + f.^q(3)); %knee 1/f model, log10 power
fitopts = optimset('Display','off','MaxFunEvals',3000,'MaxIter',3000);
apfit1 = @(lp,q0) fminsearch(@(q) sum((lp - apmodel(q,fitf)).^2), q0, fitopts); %fit every point
aprefit = @(lp,q1) fminsearch(@(q) sum((lp - apmodel(q,fitf)).^2 .* (lp - apmodel(q1,fitf) <= 0)), q1, fitopts); %refit on the points on or under the first fit q1
apstart = @(c,lp) [lp(end) + log10(8^(-c(1)) + fitf(end)^(-c(1))), sqrt(8^(-c(1))), -c(1)]; %start from a straight line (c = polyfit of log-log): knee at 8 Hz, top end lined up

%pooled fits first (straight-line start)
appq_pool = cell(nmice,npairs); %{animal x pair} [channels x 3]
for m = 1:nmice
    for j = 1:npairs
        if isempty(specpool{m,j})
            continue
        end
        appq_pool{m,j} = nan(size(specpool{m,j},1),3);
        for ch = 1:size(specpool{m,j},1)
            lp = specpool{m,j}(ch,:)'; %log10 power of this spectrum
            if all(isfinite(lp)) %flat channels stay NaN
                appq_pool{m,j}(ch,:) = aprefit(lp, apfit1(lp, apstart(polyfit(log10(fitf),lp,1),lp)));
            end
        end
    end
end

%then every trial, P2 and pre-stim (started from its channel's pooled fit, so they settle fast)
appq_tr = cell(nmice,ncond); %{animal x condition} [trials x channels x 3] P2 fits
appq_trbase = cell(nmice,ncond); %{animal x condition} [trials x channels x 3] pre-stim fits
for m = 1:nmice
    nch = numel(chanlist{m});
    for con = 1:ncond
        if ntrials(m,con) == 0 || isempty(appq_pool{m,condpair(con)})
            continue
        end
        for kind = 1:2
            if kind == 1
                S = specP{m,con};
            else
                S = specPbase{m,con};
            end
            Q = nan(ntrials(m,con),nch,3);
            for ch = 1:nch
                q0 = appq_pool{m,condpair(con)}(ch,:);
                if isnan(q0(1))
                    continue
                end
                for tr = 1:ntrials(m,con)
                    lp = double(S(tr,:,ch))';
                    if all(isfinite(lp))
                        Q(tr,ch,:) = aprefit(lp, apfit1(lp,q0));
                    end
                end
            end
            if kind == 1
                appq_tr{m,con} = Q;
            else
                appq_trbase{m,con} = Q;
            end
        end
    end
end
clear S Q

%each channel's shared band (biggest bump of the pooled spectrum above its background)
specband = cell(nmice,npairs); %{animal x pair} [channels x 1] centre of each channel's shared band, Hz
for m = 1:nmice
    for j = 1:npairs
        if isempty(appq_pool{m,j})
            continue
        end
        specband{m,j} = nan(numel(chanlist{m}),1);
        for ch = 1:numel(chanlist{m})
            if isnan(appq_pool{m,j}(ch,1)) == false
                r = specpool{m,j}(ch,:)' - apmodel(appq_pool{m,j}(ch,:),fitf); %pooled spectrum above its background
                ro = r(oscm);
                [~,im] = max(ro);
                specband{m,j}(ch) = fosc(im);
            end
        end
        if any(specband{m,j} == fosc(1) | specband{m,j} == fosc(end))
            run_notes(end+1) = mice(m) + ", " + pairnames(j) + ": " + sum(specband{m,j} == fosc(1) | specband{m,j} == fosc(end)) + ...
                " channels have their spectrum band at the edge of osc_hz (no real bump?)"; %#ok<SAGROW>
        end
    end
end

%every trial's oscillation power in its channel's shared band, P2 and pre-stim (dB above the trial's own background)
oscP2 = cell(nmice,ncond); %{animal x condition} [trials x channels] dB
oscBase = cell(nmice,ncond); %{animal x condition} [trials x channels] dB
for m = 1:nmice
    nch = numel(chanlist{m});
    for con = 1:ncond
        if isempty(appq_tr{m,con})
            continue
        end
        j = condpair(con);
        oscP2{m,con} = nan(ntrials(m,con),nch);
        oscBase{m,con} = nan(ntrials(m,con),nch);
        for ch = 1:nch
            inband = abs(fitf - specband{m,j}(ch)) <= band_hw; %all false (-> NaN) if the channel has no band
            for tr = 1:ntrials(m,con)
                q = reshape(appq_tr{m,con}(tr,ch,:),1,3);
                if isnan(q(1)) == false
                    oscP2{m,con}(tr,ch) = 10*mean(double(specP{m,con}(tr,inband,ch))' - apmodel(q,fitf(inband)));
                end
                q = reshape(appq_trbase{m,con}(tr,ch,:),1,3);
                if isnan(q(1)) == false
                    oscBase{m,con}(tr,ch) = 10*mean(double(specPbase{m,con}(tr,inband,ch))' - apmodel(q,fitf(inband)));
                end
            end
        end
    end
end
clear specPbase
elapsed_time = toc;
fprintf('Spectrum execution time: %.2f seconds\n', elapsed_time);


%% Data Analysis: P1 response of each trial (template match vs fake no-stim windows)

%Does the trial have a P1 response at all? Each channel's P1 (P1_ms, minus the mean of the p1_base_ms before the
%stimulus) is held up against that channel's average P1 waveform (the template), made from the animal's OTHER
%trials (all its loaded conditions pooled, so the trial can't match itself). The same is done at fake stimulus
%times in the pre-stim baseline (fake_ms): how well the template matches when nothing happened. It's the SHAPE
%that counts: a big spontaneous event only lines up with the template by chance (either sign), whereas just
%counting big deflections got fooled by probe-wide spontaneous events. Matches are in units of the no-stim
%spread (z); responsive channels = typical (median) trial z >= chan_z_min; trial score = mean z over them;
%P1 response = beats p1_null_prct of the fake windows' scores
disp('P1: matching every trial to its channels'' P1 templates')
tic
p1score = cell(nmice,ncond); %{animal x condition} [trials x 1] mean template z over the responsive channels
p1resp = cell(nmice,ncond); %{animal x condition} [trials x 1] true = P1 response
p1null = cell(1,nmice); %{animal} [fake windows x 1] the same score on the fake no-stim windows
p1thr = nan(1,nmice); %p1_null_prct of p1null
p1chans = cell(1,nmice); %{animal} [channels x 1] true = responsive channel
for m = 1:nmice
    cs = find(ntrials(m,:) > 0);
    for con = cs
        p1score{m,con} = nan(ntrials(m,con),1);
        p1resp{m,con} = false(ntrials(m,con),1);
    end
    if isempty(cs)
        continue
    end
    nch = numel(chanlist{m});
    W = double(cat(1, p1w{m,cs})); %[all the animal's trials x P1 samples x channels]
    Fk = double(cat(1, p1fake{m,cs})); %[all its fake windows x P1 samples x channels]
    Wsum = sum(W,1); %for the leave-one-out templates
    z = zeros(size(W,1),nch); %[trials x channels]
    zn = zeros(size(Fk,1),nch); %[fake windows x channels]
    for ch = 1:nch
        tm = Wsum(1,:,ch); %template from every trial, for the fake windows
        nproj = Fk(:,:,ch) * (tm / norm(tm))'; %[fake windows x 1]
        tmloo = Wsum(1,:,ch) - W(:,:,ch); %[trials x P1 samples] each trial's template without itself
        proj = sum(W(:,:,ch) .* tmloo, 2) ./ vecnorm(tmloo,2,2); %[trials x 1]
        z(:,ch) = proj / std(nproj); %in units of the no-stim spread
        zn(:,ch) = nproj / std(nproj);
    end
    p1chans{m} = median(z,1)' >= chan_z_min;
    if any(p1chans{m}) == false
        run_notes(end+1) = mice(m) + ": no channel with a typical P1 z >= chan_z_min, no P1 responses called"; %#ok<SAGROW>
        continue
    end
    sc = mean(z(:,p1chans{m}),2);
    p1null{m} = mean(zn(:,p1chans{m}),2);
    p1thr(m) = prctile(p1null{m}, p1_null_prct);
    %back to the conditions (same order they were stacked in)
    row = 0;
    for con = cs
        idx = row + (1:ntrials(m,con));
        p1score{m,con} = sc(idx);
        p1resp{m,con} = sc(idx) > p1thr(m);
        row = row + ntrials(m,con);
    end
end
clear p1w p1fake W Fk Wsum tmloo
elapsed_time = toc;
fprintf('P1 execution time: %.2f seconds\n', elapsed_time);


%% Data Analysis: trial classes (rhythmic response? P1 response?), and the rhythm filter

%Rhythmic RESPONSE: per animal x pair, each channel's P2 oscillation power (oscP2) in z units of its pre-stim
%oscillation power (oscBase: same band, same window length, all the pair's L and LW trials pooled). Oscillating
%channels = typical (median) trial z >= chan_z_min, L + LW pooled. Trial score = mean z over them. The same score
%on every pre-stim window is the no-stimulus reference: a trial is rhythmic when it beats rhythm_null_prct of those
rhythmscore = cell(nmice,ncond); %{animal x condition} [trials x 1] mean z over the oscillating channels
rhythmic = cell(nmice,ncond); %{animal x condition} [trials x 1] true = rhythmic response
rhythmnull = cell(nmice,npairs); %{animal x pair} [pre-stim windows x 1] the same score on the pre-stim windows
rhythmthr = nan(nmice,npairs); %rhythm_null_prct of rhythmnull
oscchans = cell(nmice,npairs); %{animal x pair} [channels x 1] true = oscillating channel
for m = 1:nmice
    for con = find(ntrials(m,:) > 0)
        rhythmscore{m,con} = nan(ntrials(m,con),1);
        rhythmic{m,con} = false(ntrials(m,con),1);
    end
    for j = 1:npairs
        if isempty(appq_pool{m,j}) %pair is missing L or LW
            continue
        end
        cs = [grpconds{j,:}];
        cs = cs(ntrials(m,cs) > 0);
        A = cat(1, oscP2{m,cs}); %[trials of the pair x channels] dB, P2
        B = cat(1, oscBase{m,cs}); %[trials of the pair x channels] dB, pre-stim
        mu = mean(B,1,'omitnan');
        sd = std(B,0,1,'omitnan');
        zA = (A - mu) ./ sd; %P2 in pre-stim z units
        zB = (B - mu) ./ sd;
        oscchans{m,j} = median(zA,1,'omitnan')' >= chan_z_min;
        if any(oscchans{m,j}) == false
            run_notes(end+1) = mice(m) + ", " + pairnames(j) + ": no channel whose P2 rhythm typically beats its pre-stim by chan_z_min, no rhythmic trials"; %#ok<SAGROW>
            continue
        end
        sc = mean(zA(:,oscchans{m,j}),2,'omitnan');
        rhythmnull{m,j} = mean(zB(:,oscchans{m,j}),2,'omitnan');
        rhythmthr(m,j) = prctile(rhythmnull{m,j}, rhythm_null_prct);
        row = 0;
        for con = cs
            idx = row + (1:ntrials(m,con));
            rhythmscore{m,con} = sc(idx);
            rhythmic{m,con} = sc(idx) > rhythmthr(m,j);
            row = row + ntrials(m,con);
        end
    end
end
clear A B zA zB

%one row per good trial (all of them, before the filter): which class is it?
classnames = ["Rhythmic+P1","Rhythmic only","P1 only","Neither"];
tc_an = strings(0,1); tc_cn = strings(0,1); tc_lvslw = strings(0,1); tc_pair = strings(0,1);
tc_fnum = []; tc_tr = []; tc_p1s = []; tc_p1 = []; tc_rs = []; tc_rh = [];
for m = 1:nmice
    for con = find(ntrials(m,:) > 0)
        n = ntrials(m,con);
        tc_an = [tc_an; repmat(mice(m),n,1)]; %#ok<AGROW>
        tc_cn = [tc_cn; repmat(string(condis{con}),n,1)]; %#ok<AGROW>
        if condgrp(con) == 1
            tc_lvslw = [tc_lvslw; repmat(condi_L,n,1)]; %#ok<AGROW>
        else
            tc_lvslw = [tc_lvslw; repmat(condi_LW,n,1)]; %#ok<AGROW>
        end
        tc_pair = [tc_pair; repmat(pairnames(condpair(con)),n,1)]; %#ok<AGROW>
        tc_fnum = [tc_fnum; trinfo{m,con}(:,1)]; %#ok<AGROW>
        tc_tr = [tc_tr; trinfo{m,con}(:,2)]; %#ok<AGROW>
        tc_p1s = [tc_p1s; p1score{m,con}]; %#ok<AGROW>
        tc_p1 = [tc_p1; p1resp{m,con}]; %#ok<AGROW>
        tc_rs = [tc_rs; rhythmscore{m,con}]; %#ok<AGROW>
        tc_rh = [tc_rh; rhythmic{m,con}]; %#ok<AGROW>
    end
end
tc_class = classnames(1 + 2*(tc_rh == 0) + (tc_p1 == 0))'; %rhythmic+P1, rhythmic only, P1 only, neither
trialclass = table(tc_an,tc_cn,tc_fnum,tc_tr,tc_lvslw,tc_pair,tc_p1s,double(tc_p1),tc_rs,double(tc_rh),tc_class, ...
    'VariableNames',["Animal_Name","Condition_Name","Condition_FileNum","Trial","LvsLW","Pair","P1_Score","P1_Resp", ...
    "Rhythm_Score","Rhythmic","TrialClass"]);
clear tc_an tc_cn tc_lvslw tc_pair tc_fnum tc_tr tc_p1s tc_p1 tc_rs tc_rh tc_class

%L vs LW: how often rhythmic, out of all trials and out of the P1-response trials only (the second one doesn't
%care how often the animal responds at all). Fisher exact per animal x pair
sk_pair = strings(0,1); sk_animal = strings(0,1); sk_counts = zeros(0,10);
sk_pRall = []; sk_pRp1 = [];
for j = 1:npairs
    for m = 1:nmice
        rows = trialclass.Animal_Name == mice(m) & trialclass.Pair == pairnames(j);
        isL = trialclass.LvsLW(rows) == condi_L;
        rh = trialclass.Rhythmic(rows) == 1;
        p1 = trialclass.P1_Resp(rows) == 1;
        if any(isL) == false || all(isL)
            continue
        end
        %counts: n, rhythmic, P1, both, neither, for L then LW
        cL = [sum(isL) sum(rh & isL) sum(p1 & isL) sum(rh & p1 & isL) sum(~rh & ~p1 & isL)];
        cLW = [sum(~isL) sum(rh & ~isL) sum(p1 & ~isL) sum(rh & p1 & ~isL) sum(~rh & ~p1 & ~isL)];
        sk_pair(end+1,1) = pairnames(j); %#ok<SAGROW>
        sk_animal(end+1,1) = mice(m); %#ok<SAGROW>
        sk_counts(end+1,:) = [cL cLW]; %#ok<SAGROW>
        [~,sk_pRall(end+1,1)] = fishertest([cL(2) cL(1)-cL(2); cLW(2) cLW(1)-cLW(2)]); %#ok<SAGROW> rhythmic out of all trials
        sk_pRp1(end+1,1) = NaN; %#ok<SAGROW>
        if cL(3) > 0 && cLW(3) > 0
            [~,sk_pRp1(end)] = fishertest([cL(4) cL(3)-cL(4); cLW(4) cLW(3)-cLW(4)]); %rhythmic out of the P1-response trials
        end
    end
end
stats_class = table(sk_pair,sk_animal,sk_counts(:,1),sk_counts(:,6),sk_counts(:,2),sk_counts(:,7),sk_counts(:,3),sk_counts(:,8), ...
    sk_counts(:,4),sk_counts(:,9),sk_counts(:,5),sk_counts(:,10), ...
    round(100*sk_counts(:,2)./sk_counts(:,1),1),round(100*sk_counts(:,7)./sk_counts(:,6),1),sk_pRall, ...
    round(100*sk_counts(:,4)./sk_counts(:,3),1),round(100*sk_counts(:,9)./sk_counts(:,8),1),sk_pRp1, ...
    'VariableNames',["Pair","Animal_Name","n_L","n_LW","nRhythmic_L","nRhythmic_LW","nP1_L","nP1_LW","nBoth_L","nBoth_LW", ...
    "nNeither_L","nNeither_LW","pctRhythmic_L","pctRhythmic_LW","p_fisher_rhythmic","pctRhythmicOfP1_L","pctRhythmicOfP1_LW", ...
    "p_fisher_rhythmicOfP1"]);
clear sk_pair sk_animal sk_counts sk_pRall sk_pRp1

%across animals (>= 2): rhythmic or not on the P1-response trials, logistic mixed model with animal as random effect
gl_pair = strings(0,1); gl_nan = []; gl_ntr = []; gl_or = []; gl_lo = []; gl_hi = []; gl_p = [];
for j = 1:npairs
    T = trialclass(trialclass.Pair == pairnames(j) & trialclass.P1_Resp == 1, :);
    if numel(unique(T.Animal_Name)) < 2 || numel(unique(T.Rhythmic)) < 2 || numel(unique(T.LvsLW)) < 2
        continue
    end
    T.Cond = categorical(T.LvsLW,[condi_L condi_LW]); %L is the reference
    T.Animal = categorical(T.Animal_Name);
    glm = fitglme(T,'Rhythmic ~ Cond + (1|Animal)','Distribution','Binomial');
    [~,~,fe] = fixedEffects(glm);
    row = strcmp(fe.Name,['Cond_' char(condi_LW)]);
    gl_pair(end+1,1) = pairnames(j); %#ok<SAGROW>
    gl_nan(end+1,1) = numel(unique(T.Animal_Name)); %#ok<SAGROW>
    gl_ntr(end+1,1) = height(T); %#ok<SAGROW>
    gl_or(end+1,1) = exp(fe.Estimate(row)); %#ok<SAGROW> odds ratio LW / L
    gl_lo(end+1,1) = exp(fe.Lower(row)); %#ok<SAGROW>
    gl_hi(end+1,1) = exp(fe.Upper(row)); %#ok<SAGROW>
    gl_p(end+1,1) = fe.pValue(row); %#ok<SAGROW>
end
stats_class_glme = table(gl_pair,gl_nan,gl_ntr,gl_or,gl_lo,gl_hi,gl_p, ...
    'VariableNames',["Pair","n_animals","n_P1trials","OddsRatio_LWvsL","OR_CIlo","OR_CIhi","p"]);
clear gl_pair gl_nan gl_ntr gl_or gl_lo gl_hi gl_p T

%the filter: keep only the rhythmic trials (rhythm_filter = 1). an_rh / an_p1 = the classes of the analyzed trials (for the table)
ntrials_all = ntrials; %[animals x conditions] good trials before the filter
an_rh = rhythmic; %{animal x condition} [analyzed trials x 1]
an_p1 = p1resp;
if rhythm_filter == 1
    for m = 1:nmice
        for con = find(ntrials_all(m,:) > 0)
            keep = rhythmic{m,con};
            lfp{m,con} = lfp{m,con}(keep,:,:);
            trinfo{m,con} = trinfo{m,con}(keep,:);
            specP{m,con} = specP{m,con}(keep,:,:);
            if isempty(appq_tr{m,con}) == false
                appq_tr{m,con} = appq_tr{m,con}(keep,:,:);
                oscP2{m,con} = oscP2{m,con}(keep,:);
            end
            an_rh{m,con} = rhythmic{m,con}(keep);
            an_p1{m,con} = p1resp{m,con}(keep);
            ntrials(m,con) = sum(keep);
            if ntrials(m,con) < 2
                run_notes(end+1) = mice(m) + ": only " + ntrials(m,con) + " rhythmic trials in " + string(condis{con}) + ", too few for most statistics"; %#ok<SAGROW>
            end
        end
    end
    %induced-only version again, from the analyzed trials (remove_erp = 1)
    if remove_erp == 1
        for m = 1:nmice
            for con = 1:ncond
                if ntrials(m,con) > 0
                    lfp_ind{m,con} = lfp{m,con} - mean(lfp{m,con},1);
                else
                    lfp_ind{m,con} = lfp{m,con};
                end
            end
        end
    else
        lfp_ind = lfp;
    end
end


%% Data Analysis 1: autocorrelation of each trial in the P2 window

disp('ACF: autocorrelating every trial and channel')
tic
lags_ms = 0:acf_maxlag_ms; %1 kHz, so lag in samples = ms
nlag = numel(lags_ms);
nP2 = numel(P2idx);
nfft = 2^nextpow2(2*nP2-1); %long enough that the fft trick gives the plain (not wrapped-around) autocorrelation

acf = cell(nmice,ncond); %{animal x condition} [trials x lags x channels], 1 at lag 0
for m = 1:nmice
    for con = 1:ncond
        if ntrials(m,con) == 0
            continue
        end
        nch = numel(chanlist{m});
        %acf{m,con} = zeros(ntrials(m,con),nlag,nch); %commented out: double, ~1 GB for 5 animals x 4 conditions
        acf{m,con} = zeros(ntrials(m,con),nlag,nch,'single'); %half the memory, plenty precise for a correlation
        for ch = 1:nch
            x = double(lfp_ind{m,con}(:,P2idx,ch)); %[trials x P2 samples]
            x = x - mean(x,2); %demean each trial's window
            F = fft(x,nfft,2); %along time
            ac = real(ifft(abs(F).^2,[],2)); %autocorrelation from the power spectrum (= xcorr, 'biased')
            acf{m,con}(:,:,ch) = ac(:,1:nlag) ./ ac(:,1); %lag 0 = 1 (= xcorr 'coeff')
        end
    end
end

%mean ACF trace of each channel, per pair and group (L / LW)
acfmean = cell(nmice,npairs,2); %{animal x pair x group} [lags x channels]
acfsem = cell(nmice,npairs,2);
for m = 1:nmice
    for j = 1:npairs
        for g = 1:2
            A = cat(1, acf{m,grpconds{j,g}}); %[trials of the group x lags x channels]
            if isempty(A)
                continue
            end
            acfmean{m,j,g} = reshape(mean(A,1), nlag, []);
            acfsem{m,j,g} = reshape(std(A,0,1), nlag, []) / sqrt(size(A,1));
        end
    end
end


%% Data Analysis 1b: autocorrelation metrics, per trial and on each channel's mean trace

%every trace gets the same 5 numbers, so they're stacked into blocks of [traces x lags] and one loop does them all:
%a block per animal x condition (every trial x channel), and one per animal x pair x group (the channel means)
%   1 decay: first lag below 1/e (ms)
%   2 area: sum over lags 1..max (ms)
%   3 secondary peak height: first peak after lag 0 sticking out by >= pk_minprom
%   4 secondary peak lag (ms)
%   5 damping: minus the slope of log(height) of the positive peaks vs lag, per second
%metricnames = ["P2_ACFDecay","P2_ACFArea","P2_ACFPk2","P2_ACFPk2Lag","P2_ACFPkDamp","All_TrialCorr","All_ChanCorr"]; %commented out: spectrum metrics (section 1c) added on the end
metricnames = ["P2_ACFDecay","P2_ACFArea","P2_ACFPk2","P2_ACFPk2Lag","P2_ACFPkDamp","All_TrialCorr","All_ChanCorr", ...
    "P2_SpecExp","P2_SpecOscPow","P2_SpecPkFreq"];
nacfmet = 5;
specidx = 8:10; %where the spectrum metrics sit in metricnames

blocks = {};
blockinfo = zeros(0,4); %[kind m con-or-pair group], kind 1 = trials of a condition, 2 = channel means of a pair group
for m = 1:nmice
    for con = 1:ncond
        if ntrials(m,con) > 0
            blocks{end+1} = reshape(permute(acf{m,con},[1 3 2]), [], nlag); %#ok<SAGROW> rows = trial + (channel-1)*trials
            blockinfo(end+1,:) = [1 m con 0]; %#ok<SAGROW>
        end
    end
    for j = 1:npairs
        for g = 1:2
            if isempty(acfmean{m,j,g}) == false
                blocks{end+1} = acfmean{m,j,g}'; %#ok<SAGROW> rows = channels
                blockinfo(end+1,:) = [2 m j g]; %#ok<SAGROW>
            end
        end
    end
end

acfmet = cell(nmice,ncond); %{animal x condition} [trials x channels x 5]
acfmet_mean = cell(nmice,npairs,2); %{animal x pair x group} [channels x 5]
for b = 1:numel(blocks)
    Y = blocks{b};
    met = nan(size(Y,1),nacfmet);
    for r = 1:size(Y,1)
        y = Y(r,:);
        if any(isnan(y)) %flat channel (0/0)
            continue
        end
        i1e = find(y < exp(-1),1);
        if isempty(i1e) == false
            met(r,1) = lags_ms(i1e); %decay, ms
        end
        met(r,2) = sum(y(2:end)); %area, ms (1 ms per lag)
        %[pks,locs] = findpeaks(y,'MinPeakProminence',pk_minprom); %commented out: chronux has its own findpeaks that shadows the Signal Processing one
        ispk = islocalmax(y,'MinProminence',pk_minprom); %peaks that stick out by at least pk_minprom
        ispk(1) = false; %lag 0 is never a secondary peak
        locs = find(ispk);
        pks = y(locs);
        if isempty(pks) == false
            met(r,3) = pks(1); %first secondary peak
            met(r,4) = lags_ms(locs(1));
            pos = pks > 0; %log needs positive heights
            if sum(pos) >= 2
                c = polyfit(lags_ms(locs(pos)), log(pks(pos)), 1); %log(height) = c(1)*lag + c(2)
                met(r,5) = -c(1) * 1000; %per second, bigger = dies out faster
            end
        end
    end

    m = blockinfo(b,2);
    if blockinfo(b,1) == 1
        con = blockinfo(b,3);
        acfmet{m,con} = reshape(met, ntrials(m,con), [], nacfmet); %[trials x channels x metrics]
    else
        acfmet_mean{m,blockinfo(b,3),blockinfo(b,4)} = met; %[channels x metrics]
    end
end
clear blocks Y
elapsed_time = toc;
fprintf('ACF execution time: %.2f seconds\n', elapsed_time);


%% Data Analysis 1c, part 2: spectrum numbers of the analyzed trials (rhythmic only when rhythm_filter = 1)

%part 1 (before the trial classes) did every trial's spectrum, its background fit and the shared band. Here: each
%group's mean spectrum (mean in power, then log) of the analyzed trials and its background fit (started from the
%channel's pooled fit), and the per-trial numbers of the analyzed trials
disp('Spectrum: group spectra of the analyzed trials')
specmeanP = cell(nmice,npairs,2); %{animal x pair x group} [channels x fitf] log10 power
appq_mean = cell(nmice,npairs,2); %{animal x pair x group} [channels x 3]
for m = 1:nmice
    nch = numel(chanlist{m});
    for j = 1:npairs
        if isempty(appq_pool{m,j})
            continue
        end
        for g = 1:2
            S = cat(1, specP{m,grpconds{j,g}}); %[analyzed trials of the group x fitf x channels]
            if isempty(S)
                continue
            end
            specmeanP{m,j,g} = log10(reshape(mean(10.^double(S),1),nfit,[]))';
            appq_mean{m,j,g} = nan(nch,3);
            for ch = 1:nch
                lp = specmeanP{m,j,g}(ch,:)';
                q0 = appq_pool{m,j}(ch,:);
                if all(isfinite(lp)) && isnan(q0(1)) == false
                    appq_mean{m,j,g}(ch,:) = aprefit(lp, apfit1(lp,q0));
                end
            end
        end
    end
end
clear S

%the numbers: exponent, oscillation power in the shared band (dB), peak frequency
specmet_mean = cell(nmice,npairs,2); %{animal x pair x group} [channels x 3] exponent, osc power (dB), group's own peak (Hz)
specmet = cell(nmice,ncond); %{animal x condition} [trials x channels x 3] exponent, osc power (dB), shared band centre (Hz)
specres = cell(nmice,ncond); %{animal x condition} [trials x fitf x channels] dB above the trial's own background (plots)
for m = 1:nmice
    nch = numel(chanlist{m});
    for j = 1:npairs
        for g = 1:2
            if isempty(appq_mean{m,j,g})
                continue
            end
            specmet_mean{m,j,g} = nan(nch,3);
            for ch = 1:nch
                q = appq_mean{m,j,g}(ch,:);
                if isnan(q(1)) || isnan(specband{m,j}(ch))
                    continue
                end
                r = specmeanP{m,j,g}(ch,:)' - apmodel(q,fitf);
                inband = abs(fitf - specband{m,j}(ch)) <= band_hw;
                ro = r(oscm);
                [~,im] = max(ro);
                specmet_mean{m,j,g}(ch,:) = [q(3) 10*mean(r(inband)) fosc(im)];
            end
        end
    end
    for con = 1:ncond
        j = condpair(con);
        if isempty(appq_tr{m,con}) || ntrials(m,con) == 0
            continue
        end
        ntr = ntrials(m,con);
        specmet{m,con} = nan(ntr,nch,3);
        specres{m,con} = nan(ntr,nfit,nch,'single');
        for ch = 1:nch
            for tr = 1:ntr
                q = reshape(appq_tr{m,con}(tr,ch,:),1,3);
                if isnan(q(1))
                    continue
                end
                specres{m,con}(tr,:,ch) = 10*(double(specP{m,con}(tr,:,ch))' - apmodel(q,fitf));
                specmet{m,con}(tr,ch,:) = [q(3) oscP2{m,con}(tr,ch) specband{m,j}(ch)];
            end
        end
    end
end


%% Data Analysis 2: trial-to-trial correlation within each channel (corr_ms window)

%For each animal x pair x channel, every trial is correlated with every other trial (L and LW together,
%one corr call). Mean pairwise r of a group = Fisher z mean of its trial pairs, back to r (the mean r
%doesn't depend on trial count, unlike ITPC). Test: shuffle the L/LW labels n_perm times on the same
%matrix and see how often the shuffled LW - L is as big as the real one (two-sided)
disp('TrialCorr: correlating trials with each other')
tic
rng('default') %same shuffles every run

trialcorr = cell(nmice,ncond); %{animal x condition} [trials x channels]
tcmean = cell(nmice,npairs); %{animal x pair} [channels x 2] mean pairwise r, L and LW
tc_p = cell(nmice,npairs); %{animal x pair} [channels x 1] permutation p
for m = 1:nmice
    nch = numel(chanlist{m});
    for con = 1:ncond
        trialcorr{m,con} = nan(ntrials(m,con),nch);
    end
    for j = 1:npairs
        tcmean{m,j} = nan(nch,2);
        tc_p{m,j} = nan(nch,1);
        cs = [grpconds{j,1} grpconds{j,2}];
        if sum(ntrials(m,grpconds{j,1})) < 2 || sum(ntrials(m,grpconds{j,2})) < 2
            run_notes(end+1) = mice(m) + ", " + pairnames(j) + ": fewer than 2 trials in L or LW, no trial-trial correlation"; %#ok<SAGROW>
            continue
        end

        %which condition and trial each row of the stacked trials is
        rowcond = [];
        rowtr = [];
        for c = cs
            rowcond = [rowcond; repmat(c,ntrials(m,c),1)]; %#ok<AGROW>
            rowtr = [rowtr; (1:ntrials(m,c))']; %#ok<AGROW>
        end
        isLW = condgrp(rowcond)' == 2;
        ntot = numel(rowcond);
        X = cat(1, lfp{m,cs}); %[all trials of the pair x keep_ms x channels], NOT erp-removed

        for ch = 1:nch
            C = corr(double(X(:,corridx,ch))'); %[trials x trials]
            C(C > 0.9999) = 0.9999; %keep atanh finite
            Z = atanh(C); %Fisher z
            Z(1:ntot+1:end) = NaN; %no trial with itself

            %each trial's mean z with the other trials of its own group
            zrowL = mean(Z(:,isLW == false),2,'omitnan');
            zrowLW = mean(Z(:,isLW),2,'omitnan');
            own = zrowL;
            own(isLW) = zrowLW(isLW);
            for c = cs
                trialcorr{m,c}(rowtr(rowcond == c),ch) = tanh(own(rowcond == c));
            end

            zL = mean(Z(isLW == false,isLW == false),'all','omitnan'); %mean over L trial pairs
            zLW = mean(Z(isLW,isLW),'all','omitnan');
            tcmean{m,j}(ch,:) = tanh([zL zLW]);

            obs = zLW - zL;
            null = zeros(n_perm,1);
            for k = 1:n_perm
                lab = isLW(randperm(ntot)); %shuffled labels
                null(k) = mean(Z(lab,lab),'all','omitnan') - mean(Z(lab == false,lab == false),'all','omitnan');
            end
            if isnan(obs) == false
                tc_p{m,j}(ch) = (sum(abs(null) >= abs(obs)) + 1) / (n_perm + 1);
            end
        end
    end
end
clear X C Z
elapsed_time = toc;
fprintf('TrialCorr execution time: %.2f seconds\n', elapsed_time);


%% Data Analysis 3: channel-to-channel correlation on each trial (corr_ms window)

disp('ChanCorr: correlating channels with each other')
tic
chcorr = cell(nmice,ncond); %{animal x condition} [trials x channels] channel's mean r with the other channels
offdiag = cell(nmice,ncond); %{animal x condition} [trials x 1] mean r over all channel pairs
chanZsum = cell(nmice,ncond); %{animal x condition} [channels x channels] sum of Fisher z over trials
for m = 1:nmice
    nch = numel(chanlist{m});
    upmask = triu(true(nch),1); %each channel pair once
    for con = 1:ncond
        ntr = ntrials(m,con);
        chcorr{m,con} = nan(ntr,nch);
        offdiag{m,con} = nan(ntr,1);
        chanZsum{m,con} = zeros(nch);
        for tr = 1:ntr
            C = corr(reshape(double(lfp_ind{m,con}(tr,corridx,:)), [], nch)); %[channels x channels]
            C(C > 0.9999) = 0.9999; %keep atanh finite
            Z = atanh(C);
            Z(1:nch+1:end) = NaN; %no channel with itself
            chcorr{m,con}(tr,:) = tanh(mean(Z,2,'omitnan'))';
            offdiag{m,con}(tr) = tanh(mean(Z(upmask),'omitnan'));
            chanZsum{m,con} = chanZsum{m,con} + Z;
        end
    end
end
clear C Z
elapsed_time = toc;
fprintf('ChanCorr execution time: %.2f seconds\n', elapsed_time);


%% Data PreProcessing: one table, one row per trial x channel

%rows go animal -> condition -> channel -> trial (trial fastest), like the summary csv
col_an = strings(0,1); col_cn = strings(0,1); col_reg = strings(0,1); col_lvslw = strings(0,1); col_pair = strings(0,1);
col_fnum = []; col_tr = []; col_ch = []; col_shk = []; col_dep = []; col_rh = []; col_p1 = [];
col_rank = zeros(0,numel(shankscore_use));
col_met = zeros(0,numel(metricnames));
for m = 1:nmice
    nch = numel(chanlist{m});
    for con = 1:ncond
        ntr = ntrials(m,con);
        if ntr == 0
            continue
        end
        nrow = ntr*nch;
        col_an = [col_an; repmat(mice(m),nrow,1)]; %#ok<AGROW>
        col_cn = [col_cn; repmat(string(condis{con}),nrow,1)]; %#ok<AGROW>
        col_reg = [col_reg; repmat(region,nrow,1)]; %#ok<AGROW>
        if condgrp(con) == 1
            col_lvslw = [col_lvslw; repmat(condi_L,nrow,1)]; %#ok<AGROW>
        else
            col_lvslw = [col_lvslw; repmat(condi_LW,nrow,1)]; %#ok<AGROW>
        end
        col_pair = [col_pair; repmat(pairnames(condpair(con)),nrow,1)]; %#ok<AGROW>
        col_fnum = [col_fnum; repmat(trinfo{m,con}(:,1),nch,1)]; %#ok<AGROW>
        col_tr = [col_tr; repmat(trinfo{m,con}(:,2),nch,1)]; %#ok<AGROW>
        col_ch = [col_ch; repelem(chanlist{m},ntr,1)]; %#ok<AGROW> each channel ntr times in a row
        col_shk = [col_shk; repelem(chshank{m},ntr,1)]; %#ok<AGROW>
        col_dep = [col_dep; repelem(chdepth{m},ntr,1)]; %#ok<AGROW>
        col_rank = [col_rank; repelem(chrank{m},ntr,1)]; %#ok<AGROW>
        col_rh = [col_rh; repmat(double(an_rh{m,con}),nch,1)]; %#ok<AGROW> trial class, same for every channel of the trial
        col_p1 = [col_p1; repmat(double(an_p1{m,con}),nch,1)]; %#ok<AGROW>
        %mvals = [reshape(acfmet{m,con},nrow,nacfmet) trialcorr{m,con}(:) chcorr{m,con}(:)]; %commented out: spectrum metrics added on the end
        if isempty(specmet{m,con}) %pair is missing L or LW, no spectrum fit
            smet = nan(nrow,numel(specidx));
        else
            smet = reshape(specmet{m,con},nrow,numel(specidx));
        end
        mvals = [reshape(acfmet{m,con},nrow,nacfmet) trialcorr{m,con}(:) chcorr{m,con}(:) smet]; %[rows x metrics]
        col_met = [col_met; mvals]; %#ok<AGROW>
    end
end

tableres = table(col_an,col_cn,col_fnum,col_tr,col_reg,col_ch,col_shk,col_dep,col_lvslw,col_pair, ...
    'VariableNames',["Animal_Name","Condition_Name","Condition_FileNum","Trial","Region","Channel_ID", ...
    "Chan_Shank","Chan_Depth","LvsLW","Pair"]);
for k = 1:numel(rankcols)
    tableres.(rankcols(k)) = col_rank(:,k);
end
tableres.Rhythmic = col_rh; %1 = rhythmic response (all 1 when rhythm_filter = 1)
tableres.P1_Resp = col_p1; %1 = P1 response
for k = 1:numel(metricnames)
    tableres.(metricnames(k)) = col_met(:,k);
end
clear col_an col_cn col_fnum col_tr col_reg col_ch col_shk col_dep col_lvslw col_pair col_rank col_met mvals smet col_rh col_p1


%% Data Analysis: per channel L vs LW ("units" = animal x channel x pair) and per animal

nmet = numel(metricnames);
uL = zeros(0,nmet); uLW = zeros(0,nmet); uP = zeros(0,nmet); uQ = zeros(0,nmet); %channel-level L, LW, p, q
u_animal = strings(0,1); u_chan = []; u_pair = strings(0,1);
u_ntrL = []; u_ntrLW = []; u_pkfracL = []; u_pkfracLW = [];

%per animal (one row per pair x animal x metric)
sa_pair = strings(0,1); sa_animal = strings(0,1); sa_metric = strings(0,1);
sa_n = []; sa_medL = []; sa_medLW = []; sa_higher = []; sa_p = []; sa_nsig = [];
%channel-channel per animal (one row per pair x animal)
sc_pair = strings(0,1); sc_animal = strings(0,1);
sc_rL = []; sc_rLW = []; sc_ntrL = []; sc_ntrLW = []; sc_p = [];

for m = 1:nmice
    nch = numel(chanlist{m});
    for j = 1:npairs
        cL = grpconds{j,1};
        cLW = grpconds{j,2};
        if isempty(acfmean{m,j,1}) || isempty(acfmean{m,j,2})
            run_notes(end+1) = mice(m) + ", " + pairnames(j) + ": missing L or LW, left out of the comparisons"; %#ok<SAGROW>
            continue
        end
        AL = cat(1, acfmet{m,cL}); %[L trials x channels x 5]
        ALW = cat(1, acfmet{m,cLW});
        ChL = cat(1, chcorr{m,cL}); %[L trials x channels]
        ChLW = cat(1, chcorr{m,cLW});

        bL = nan(nch,nmet); bLW = nan(nch,nmet); bP = nan(nch,nmet); bQ = nan(nch,nmet);
        %ACF metrics: channel value from the mean ACF trace, test on the per-trial values
        for k = 1:nacfmet
            bL(:,k) = acfmet_mean{m,j,1}(:,k);
            bLW(:,k) = acfmet_mean{m,j,2}(:,k);
            for ch = 1:nch
                a = AL(:,ch,k);
                b = ALW(:,ch,k);
                a = a(isnan(a) == false);
                b = b(isnan(b) == false);
                if numel(a) >= 2 && numel(b) >= 2
                    bP(ch,k) = ranksum(a,b);
                end
            end
        end
        %spectrum: channel value from the group's mean spectrum, test on the per-trial values
        SpL = cat(1, specmet{m,cL}); %[L trials x channels x 3]
        SpLW = cat(1, specmet{m,cLW});
        for k = 1:numel(specidx)
            bL(:,specidx(k)) = specmet_mean{m,j,1}(:,k);
            bLW(:,specidx(k)) = specmet_mean{m,j,2}(:,k);
            if metricnames(specidx(k)) == "P2_SpecPkFreq" %every trial carries the shared band here, so no per-trial test
                continue
            end
            for ch = 1:nch
                a = SpL(:,ch,k);
                b = SpLW(:,ch,k);
                a = a(isnan(a) == false);
                b = b(isnan(b) == false);
                if numel(a) >= 2 && numel(b) >= 2
                    bP(ch,specidx(k)) = ranksum(a,b);
                end
            end
        end
        %trial-trial: mean pairwise r and the permutation p
        bL(:,6) = tcmean{m,j}(:,1);
        bLW(:,6) = tcmean{m,j}(:,2);
        bP(:,6) = tc_p{m,j};
        %channel-channel: mean over trials (in z), test on the per-trial values
        bL(:,7) = tanh(mean(atanh(ChL),1,'omitnan'))';
        bLW(:,7) = tanh(mean(atanh(ChLW),1,'omitnan'))';
        for ch = 1:nch
            a = ChL(isnan(ChL(:,ch)) == false,ch);
            b = ChLW(isnan(ChLW(:,ch)) == false,ch);
            if numel(a) >= 2 && numel(b) >= 2
                bP(ch,7) = ranksum(a,b);
            end
        end

        %FDR (Benjamini-Hochberg) over this animal's channels, per metric
        for k = 1:nmet
            pch = bP(:,k);
            hasp = find(isnan(pch) == false);
            [ps,ord] = sort(pch(hasp));
            nt = numel(ps);
            qs = ps .* nt ./ (1:nt)';
            qs = flipud(cummin(flipud(qs))); %make them non-decreasing from the top
            bQ(hasp(ord),k) = min(qs,1);
        end

        %stack into the unit arrays
        uL = [uL; bL]; uLW = [uLW; bLW]; uP = [uP; bP]; uQ = [uQ; bQ]; %#ok<AGROW>
        u_animal = [u_animal; repmat(mice(m),nch,1)]; %#ok<AGROW>
        u_chan = [u_chan; chanlist{m}]; %#ok<AGROW>
        u_pair = [u_pair; repmat(pairnames(j),nch,1)]; %#ok<AGROW>
        u_ntrL = [u_ntrL; repmat(size(AL,1),nch,1)]; %#ok<AGROW>
        u_ntrLW = [u_ntrLW; repmat(size(ALW,1),nch,1)]; %#ok<AGROW>
        u_pkfracL = [u_pkfracL; mean(isnan(AL(:,:,3)) == false,1)']; %#ok<AGROW> fraction of trials with a secondary peak
        u_pkfracLW = [u_pkfracLW; mean(isnan(ALW(:,:,3)) == false,1)']; %#ok<AGROW>

        %per animal: signrank across its channels, per metric
        for k = 1:nmet
            use = isnan(bL(:,k)) == false & isnan(bLW(:,k)) == false;
            sa_pair(end+1,1) = pairnames(j); %#ok<SAGROW>
            sa_animal(end+1,1) = mice(m); %#ok<SAGROW>
            sa_metric(end+1,1) = metricnames(k); %#ok<SAGROW>
            sa_n(end+1,1) = sum(use); %#ok<SAGROW>
            sa_medL(end+1,1) = median(bL(use,k)); %#ok<SAGROW>
            sa_medLW(end+1,1) = median(bLW(use,k)); %#ok<SAGROW>
            sa_higher(end+1,1) = sum(bLW(use,k) > bL(use,k)); %#ok<SAGROW>
            sa_p(end+1,1) = NaN; %#ok<SAGROW>
            if sum(use) > 0
                sa_p(end) = signrank(bLW(use,k),bL(use,k));
            end
            sa_nsig(end+1,1) = sum(bQ(:,k) < 0.05); %#ok<SAGROW>
        end

        %channel-channel: mean over all channel pairs, L vs LW over trials
        oL = cat(1, offdiag{m,cL});
        oLW = cat(1, offdiag{m,cLW});
        sc_pair(end+1,1) = pairnames(j); %#ok<SAGROW>
        sc_animal(end+1,1) = mice(m); %#ok<SAGROW>
        sc_rL(end+1,1) = tanh(mean(atanh(oL),'omitnan')); %#ok<SAGROW>
        sc_rLW(end+1,1) = tanh(mean(atanh(oLW),'omitnan')); %#ok<SAGROW>
        sc_ntrL(end+1,1) = numel(oL); %#ok<SAGROW>
        sc_ntrLW(end+1,1) = numel(oLW); %#ok<SAGROW>
        sc_p(end+1,1) = ranksum(oL,oLW); %#ok<SAGROW>
    end
end

%unit table: one row per animal x channel x pair, columns <metric>_L/_LW/_p/_q
unitres = table(u_animal,u_chan,u_pair,u_ntrL,u_ntrLW,u_pkfracL,u_pkfracLW, ...
    'VariableNames',["Animal_Name","Channel_ID","Pair","n_L","n_LW","Pk2Frac_L","Pk2Frac_LW"]);
for k = 1:nmet
    unitres.(metricnames(k) + "_L") = uL(:,k);
    unitres.(metricnames(k) + "_LW") = uLW(:,k);
    unitres.(metricnames(k) + "_p") = uP(:,k);
    unitres.(metricnames(k) + "_q") = uQ(:,k);
end
stats_animal = table(sa_pair,sa_animal,sa_metric,sa_n,sa_medL,sa_medLW,sa_higher,sa_p,sa_nsig, ...
    'VariableNames',["Pair","Animal_Name","Metric","n_channels","median_L","median_LW","n_channels_LWhigher", ...
    "p_signrank","n_channels_q05"]);
stats_chanchan = table(sc_pair,sc_animal,sc_rL,sc_rLW,sc_ntrL,sc_ntrLW,sc_p, ...
    'VariableNames',["Pair","Animal_Name","meanChanPair_r_L","meanChanPair_r_LW","n_L","n_LW","p_ranksum"]);


%% Data Analysis: across animals (paired units, per pair and all together, + mixed model)

%same idea as the P1 script: signrank across units (rough, channels of one mouse aren't independent) and a
%mixed model with each animal getting its own LW effect (1 + Cond|Animal) and channel as a random effect.
%scale per metric: log for the positive times (effect = LW/L ratio), Fisher z for correlations
%(effect = difference in z), raw for the rest (effect = LW - L). All p values two-sided
%metric_scale = ["log","raw","raw","log","raw","z","z"]; %commented out: spectrum metrics added on the end
metric_scale = ["log","raw","raw","log","raw","z","z","raw","raw","log"]; %same order as metricnames (exponent raw, dB raw, Hz log)
pairlist = [pairnames(:); "all"];
npl = numel(pairlist);
nrows = nmet*npl;

st_metric = strings(nrows,1); st_pair = strings(nrows,1); st_scale = strings(nrows,1);
n_units = zeros(nrows,1); median_L = nan(nrows,1); median_LW = nan(nrows,1);
n_units_LWhigher = zeros(nrows,1); p_signrank = nan(nrows,1);
lme_n_animals = zeros(nrows,1); lme_effect = nan(nrows,1); lme_CIlo = nan(nrows,1); lme_CIhi = nan(nrows,1); lme_p = nan(nrows,1);
lme_models = cell(nrows,1);
for k = 1:nmet
    for j = 1:npl
        r = (k-1)*npl + j;
        st_metric(r) = metricnames(k);
        st_pair(r) = pairlist(j);
        st_scale(r) = metric_scale(k);
        vL = uL(:,k);
        vLW = uLW(:,k);
        use = isnan(vL) == false & isnan(vLW) == false;
        if pairlist(j) ~= "all"
            use = use & u_pair == pairlist(j);
        end
        if metric_scale(k) == "log"
            use = use & vL > 0 & vLW > 0; %log needs positive values
        end
        n_units(r) = sum(use);
        if n_units(r) == 0
            continue
        end
        median_L(r) = median(vL(use));
        median_LW(r) = median(vLW(use));
        n_units_LWhigher(r) = sum(vLW(use) > vL(use));
        p_signrank(r) = signrank(vLW(use),vL(use));

        %mixed model, long format: one row per unit x condition
        if metric_scale(k) == "log"
            tL = log(vL(use)); tLW = log(vLW(use));
        elseif metric_scale(k) == "z"
            tL = atanh(vL(use)); tLW = atanh(vLW(use));
        else
            tL = vL(use); tLW = vLW(use);
        end
        nu = sum(use);
        lmetab = table([tL; tLW], ...
            categorical([repmat(condi_L,nu,1); repmat(condi_LW,nu,1)],[condi_L condi_LW]), ... %L is the reference
            categorical([u_pair(use); u_pair(use)]), ...
            [u_animal(use); u_animal(use)], ...
            [u_animal(use) + "_ch" + string(u_chan(use)); u_animal(use) + "_ch" + string(u_chan(use))], ...
            'VariableNames',["Value","Cond","Pair","Animal","ChanKey"]);
        lme_n_animals(r) = numel(unique(lmetab.Animal));
        if lme_n_animals(r) < 2
            if k == 1 %same for every metric, so just say it once
                run_notes(end+1) = pairlist(j) + ": only 1 animal, no mixed model"; %#ok<SAGROW>
            end
            continue
        end
        if pairlist(j) == "all" && numel(unique(lmetab.Pair)) > 1
            formula = 'Value ~ Cond + Pair + (1 + Cond|Animal) + (1|ChanKey)'; %pair as a fixed effect, LW effect averaged over pairs
        else
            formula = 'Value ~ Cond + (1 + Cond|Animal) + (1|ChanKey)';
        end
        mdl = fitlme(lmetab,formula,'FitMethod','REML'); %REML: better variance estimates with few animals
        lme_models{r} = mdl;
        [~,~,fe] = fixedEffects(mdl,'DFMethod','satterthwaite'); %DF from the animals, not the channel count
        row = strcmp(fe.Name,['Cond_' char(condi_LW)]);
        if metric_scale(k) == "log"
            lme_effect(r) = exp(fe.Estimate(row)); %LW / L ratio
            lme_CIlo(r) = exp(fe.Lower(row));
            lme_CIhi(r) = exp(fe.Upper(row));
        else
            lme_effect(r) = fe.Estimate(row); %LW - L (in z for the correlations)
            lme_CIlo(r) = fe.Lower(row);
            lme_CIhi(r) = fe.Upper(row);
        end
        lme_p(r) = fe.pValue(row);
    end
end
stats_paired = table(st_metric,st_pair,n_units,median_L,median_LW,n_units_LWhigher,p_signrank, ...
    lme_n_animals,st_scale,lme_effect,lme_CIlo,lme_CIhi,lme_p, ...
    'VariableNames',["Metric","Pair","n_units","median_L","median_LW","n_units_LWhigher","p_signrank", ...
    "lme_n_animals","lme_scale","lme_effect","lme_CIlo","lme_CIhi","lme_p"]);


%% Data Visualization: mean autocorrelation and spectrum per channel in the probe layout (L vs LW)

%one figure per animal x pair, for the ACF (plot_acf) and for the spectrum above its 1/f background
%(plot_spec, mean of each trial's dB above its own background). tiledlayout with no spacing instead of
%subplot, so each channel's panel is bigger; the channel number sits inside the panel instead of a title.
%Dots: ACF = the mean trace's secondary peak, spectrum = the channel's shared band and its oscillation power
plotkinds = find([plot_acf plot_spec] == 1); %1 = ACF, 2 = spectrum
colL = [0 0 0]; %L black
colLW = [0 0.6 0]; %LW green
for pk = plotkinds
    if pk == 1
        disp('ACF: Making some pictures for you')
    else
        disp('Spectrum: Making some pictures for you')
    end
    for m = 1:nmice
        if plotsave == 1 && isfolder([save_directory '\Figures_LFPCorr\' char(mice(m))]) == false
            mkdir([save_directory '\Figures_LFPCorr\' char(mice(m))])
        end
        pmap = pinfos{m}.ProbeMaps{mouseProbe(m)}; %[rows x shanks]
        [nrow,nshank] = size(pmap);
        for j = 1:npairs
            if isempty(acfmean{m,j,1}) || isempty(acfmean{m,j,2})
                continue
            end
            %the traces, their SEM and the dots, per group (L, LW)
            trmean = cell(1,2); %[x points x channels]
            trsem = cell(1,2);
            dotx = cell(1,2); %[channels x 1]
            doty = cell(1,2);
            for g = 1:2
                if pk == 1
                    trmean{g} = acfmean{m,j,g};
                    trsem{g} = acfsem{m,j,g};
                    dotx{g} = acfmet_mean{m,j,g}(:,4);
                    doty{g} = acfmet_mean{m,j,g}(:,3);
                else
                    R = cat(1, specres{m,grpconds{j,g}}); %[trials x fitf x channels] dB above background
                    trmean{g} = reshape(mean(R,1,'omitnan'), nfit, []);
                    trsem{g} = reshape(std(R,0,1,'omitnan'), nfit, []) / sqrt(size(R,1));
                    dotx{g} = specband{m,j};
                    doty{g} = specmet_mean{m,j,g}(:,2);
                end
            end
            clear R
            %one y range for the whole probe
            ylo = min([min(trmean{1} - trsem{1},[],'all') min(trmean{2} - trsem{2},[],'all')]);
            yhi = max([max(trmean{1} + trsem{1},[],'all') max(trmean{2} + trsem{2},[],'all')]);
            if pk == 1
                xs = lags_ms;
                ylo = floor(ylo*10)/10; %round down to 0.1
                yhi = 1;
                xrange = [0 acf_maxlag_ms];
                whattxt = "P2 autocorrelation";
                xlab = 'lag (ms)';
                ylab = 'autocorrelation';
                ftag = 'ACF';
            else
                xs = fitf';
                ylo = floor(ylo);
                yhi = ceil(yhi);
                xrange = [fitf(1) fitf(end)];
                whattxt = "P2 spectrum above 1/f";
                xlab = 'frequency (Hz)';
                ylab = 'dB above the 1/f background';
                ftag = 'Spec';
            end

            figure();
            set(gcf, 'Position', get(0, 'Screensize'));
            tl = tiledlayout(nrow,nshank,'TileSpacing','none','Padding','tight');
            %title(tl,char(mice(m) + " " + region + "-Probe P2 autocorrelation: " + pairnames(j) + ... %commented out: shared with the spectrum plot now
            %    " (" + condi_L + " black, " + condi_LW + " green, mean +/- SEM)"),'Interpreter','none')
            title(tl,char(mice(m) + " " + region + "-Probe " + whattxt + ": " + pairnames(j) + ...
                " (" + condi_L + " black, " + condi_LW + " green, mean +/- SEM)"),'Interpreter','none')
            xlabel(tl,xlab)
            ylabel(tl,ylab)
            hL = []; hLW = [];
            for row = 1:nrow
                for shank = 1:nshank
                    chid = pmap(row,shank);
                    ax = nexttile(tl,(row-1)*nshank + shank);
                    k = find(chanlist{m} == chid);
                    if isempty(k)
                        axis(ax,'off')
                        if chid > 0
                            text(ax,0.5,0.5,['CH' num2str(chid) ' removed'],'Units','normalized','HorizontalAlignment','center','FontSize',7)
                        end
                        continue
                    end
                    hold(ax,'on')
                    for g = 1:2
                        mu = trmean{g}(:,k)';
                        se = trsem{g}(:,k)';
                        if g == 1
                            col = colL;
                        else
                            col = colLW;
                        end
                        fill(ax,[xs fliplr(xs)],[mu+se fliplr(mu-se)],col,'FaceAlpha',0.25,'EdgeColor','none')
                        h = plot(ax,xs,mu,'Color',col,'LineWidth',1);
                        if g == 1
                            hL = h;
                        else
                            hLW = h;
                        end
                        if isnan(doty{g}(k)) == false
                            plot(ax,dotx{g}(k),doty{g}(k),'o','Color',col,'MarkerFaceColor',col,'MarkerSize',3)
                        end
                    end
                    yline(ax,0,':')
                    hold(ax,'off')
                    xlim(ax,xrange)
                    ylim(ax,[ylo yhi])
                    text(ax,0.97,0.93,['CH' num2str(chid)],'Units','normalized','HorizontalAlignment','right', ...
                        'VerticalAlignment','top','FontSize',7)
                    if row < nrow || shank > 1 %tick labels only on the bottom-left panel
                        set(ax,'XTickLabel',[],'YTickLabel',[])
                    end
                end
            end
            if isempty(hL) == false
                lgd = legend([hL hLW],cellstr([condi_L condi_LW]),'Orientation','horizontal');
                lgd.Layout.Tile = 'south';
            end

            if plotsave == 1
                fname = sprintf([char(mice(m)) '-' char(region) '_Probe-' ftag '_' regexprep(char(pairnames(j)),'\W','_')]);
                savefig(gcf,[save_directory '\Figures_LFPCorr\' char(mice(m)) '\' fname],'compact')
            end
        end
    end
end


%% Data Visualization: per-animal channel maps of LW - L

%One figure per metric x pair, one tile per animal laid out like its probe (rows top to bottom, columns
%= shanks), like the P1 script. Colour = LW - L of the channel's value (blue = LW lower, red = LW higher),
%* = that channel's L vs LW test significant after FDR within the animal. Dark grey = kept channel with
%no value (too few trials / no peak), black = removed or no site. Tile title = per-animal signrank
if plot_chanmaps == 1
    disp('Channel maps: Making some pictures for you')
    if plotsave == 1 && isfolder([save_directory '\Figures_LFPCorr']) == false
        mkdir([save_directory '\Figures_LFPCorr'])
    end
    bwr = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
    greyc = [0.4 0.4 0.4];
    for km = 1:numel(chanmap_metrics)
        k = find(metricnames == chanmap_metrics(km));
        if isempty(k)
            run_notes(end+1) = "chanmap_metrics: " + chanmap_metrics(km) + " is not a metric name, skipped"; %#ok<SAGROW>
            continue
        end
        for j = 1:npairs
            inpair = u_pair == pairnames(j);
            pair_animals = unique(u_animal(inpair),'stable');
            if isempty(pair_animals)
                continue
            end
            dmax = max(abs(uLW(inpair,k) - uL(inpair,k))); %one symmetric colour scale for all tiles
            if isnan(dmax) || dmax == 0
                dmax = 1;
            end

            figure('Name',metricnames(k) + ": " + pairnames(j));
            set(gcf, 'Position', get(0, 'Screensize'));
            tl = tiledlayout('flow','TileSpacing','compact');
            title(tl,metricnames(k) + ", LW - L per channel: " + pairnames(j),'Interpreter','none')
            for a = 1:numel(pair_animals)
                m = find(mice == pair_animals(a));
                inA = inpair & u_animal == pair_animals(a);
                ax = nexttile;
                pmap = pinfos{m}.ProbeMaps{mouseProbe(m)};

                M = nan(size(pmap)); %NaN = no value (grey or black)
                S = false(size(pmap)); %significant channels
                G = false(size(pmap)); %kept but no value
                chans = u_chan(inA);
                dvals = uLW(inA,k) - uL(inA,k);
                qvals = uQ(inA,k);
                for ch = 1:numel(chans)
                    loc = find(pmap == chans(ch));
                    M(loc) = dvals(ch);
                    S(loc) = qvals(ch) < 0.05;
                    G(loc) = isnan(dvals(ch));
                end

                h = imagesc(ax,M);
                set(h,'AlphaData',isnan(M) == false) %NaN cells see-through...
                set(ax,'Color','k') %...onto a black background (removed / no site)
                colormap(ax,bwr)
                clim(ax,[-dmax dmax])
                hold(ax,'on')
                image(ax,repmat(reshape(greyc,1,1,3),size(G)),'AlphaData',double(G)) %grey squares only where G
                [sr,sc] = find(S);
                text(ax,sc,sr,'*','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14,'FontWeight','bold')
                hold(ax,'off')
                axis(ax,'image')
                xticks(ax,1:size(pmap,2)); yticks(ax,1:size(pmap,1));
                xlabel(ax,'shank'); ylabel(ax,'channel (top = 1)')
                sarow = find(sa_pair == pairnames(j) & sa_animal == pair_animals(a) & sa_metric == metricnames(k));
                title(ax,{char(pair_animals(a)), sprintf('LW higher %d/%d, p = %.3g', sa_higher(sarow), sa_n(sarow), sa_p(sarow))},'Interpreter','none')
            end
            cb = colorbar(ax);
            cb.Layout.Tile = 'east';
            cb.Label.String = [char(metricnames(k)) ' LW - L'];
            cb.Label.Interpreter = 'none';

            %legend for the grey and black squares (invisible patches just to make the entries)
            hold(ax,'on')
            lg1 = patch(ax,NaN,NaN,greyc,'EdgeColor','none');
            lg2 = patch(ax,NaN,NaN,'k','EdgeColor','none');
            hold(ax,'off')
            lgd = legend([lg1 lg2],{'No value (too few trials / no peak)','Removed / no site'},'Orientation','horizontal');
            lgd.Layout.Tile = 'south';
            lgd.Color = 'w'; %otherwise it copies the black axes background
            lgd.TextColor = 'k';

            if plotsave == 1
                fname = sprintf(['chanmap-' char(metricnames(k)) '_' regexprep(char(pairnames(j)),'\W','_')]);
                savefig(gcf,[save_directory '\Figures_LFPCorr\' fname],'compact')
            end
        end
    end
end


%% Data Visualization: channel x channel correlation matrices

%one figure per animal x pair: mean r in L, in LW, LW - L (channels shank by shank, top to bottom,
%white lines between shanks), and each trial's mean over all channel pairs
if plot_chanchan == 1
    disp('ChanCorr: Making some pictures for you')
    bwr = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
    for m = 1:nmice
        if plotsave == 1 && isfolder([save_directory '\Figures_LFPCorr\' char(mice(m))]) == false
            mkdir([save_directory '\Figures_LFPCorr\' char(mice(m))])
        end
        ord = chorder{m};
        shk = chshank{m}(ord);
        bounds = find(diff(shk)) + 0.5; %between shanks
        for j = 1:npairs
            cL = grpconds{j,1};
            cLW = grpconds{j,2};
            nL = sum(ntrials(m,cL));
            nLW = sum(ntrials(m,cLW));
            if nL == 0 || nLW == 0
                continue
            end
            %mean z over all trials of the group, back to r
            nch = numel(ord);
            zL = zeros(nch,nch); zLW = zeros(nch,nch);
            for c = cL
                zL = zL + chanZsum{m,c};
            end
            for c = cLW
                zLW = zLW + chanZsum{m,c};
            end
            rL = tanh(zL(ord,ord) / nL); %[channels x channels], shank order
            rLW = tanh(zLW(ord,ord) / nLW);
            rlim = [min([rL(:); rLW(:)]) max([rL(:); rLW(:)])];
            dlim = max(abs(rLW(:) - rL(:)));
            if isnan(dlim) || dlim == 0
                dlim = 1;
            end

            figure();
            set(gcf, 'Position', get(0, 'Screensize'));
            tl = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
            title(tl,char(mice(m) + " " + region + "-Probe channel x channel correlation (" + corr_ms(1) + "-" + corr_ms(2) + ...
                " ms): " + pairnames(j)),'Interpreter','none')
            mats = {rL, rLW, rLW - rL};
            ttls = {char(condi_L + " (n = " + nL + ")"), char(condi_LW + " (n = " + nLW + ")"), 'LW - L'};
            for i = 1:3
                ax = nexttile;
                h = imagesc(ax,mats{i});
                set(h,'AlphaData',isnan(mats{i}) == false)
                axis(ax,'image')
                if i < 3
                    colormap(ax,parula)
                    clim(ax,rlim)
                else
                    colormap(ax,bwr)
                    clim(ax,[-dlim dlim])
                end
                colorbar(ax)
                hold(ax,'on')
                for bd = bounds'
                    xline(ax,bd,'w','LineWidth',1)
                    yline(ax,bd,'w','LineWidth',1)
                end
                hold(ax,'off')
                ticks = 1:numel(ord);
                set(ax,'XTick',ticks,'XTickLabel',chanlist{m}(ord),'YTick',ticks,'YTickLabel',chanlist{m}(ord),'FontSize',5)
                title(ax,ttls{i},'FontSize',10)
            end
            %each trial's mean over all channel pairs
            ax = nexttile;
            oL = cat(1, offdiag{m,cL});
            oLW = cat(1, offdiag{m,cLW});
            hold(ax,'on')
            scatter(ax,1 + 0.15*randn(numel(oL),1),oL,12,[0 0 0],'filled','MarkerFaceAlpha',0.4)
            scatter(ax,2 + 0.15*randn(numel(oLW),1),oLW,12,[0 0.6 0],'filled','MarkerFaceAlpha',0.4)
            plot(ax,[1 2],[median(oL,'omitnan') median(oLW,'omitnan')],'r-o','LineWidth',2)
            hold(ax,'off')
            xlim(ax,[0.5 2.5]); xticks(ax,[1 2]); xticklabels(ax,[condi_L condi_LW])
            ylabel(ax,'mean r over channel pairs (per trial)')
            scrow = find(sc_pair == pairnames(j) & sc_animal == mice(m));
            if isempty(scrow) == false
                title(ax,sprintf('ranksum p = %.3g', sc_p(scrow)),'FontSize',10)
            end

            if plotsave == 1
                fname = sprintf([char(mice(m)) '-' char(region) '_Probe-ChanCorr_' regexprep(char(pairnames(j)),'\W','_')]);
                savefig(gcf,[save_directory '\Figures_LFPCorr\' char(mice(m)) '\' fname],'compact')
            end
        end
    end
end


%% Data Visualization: across animals, paired L vs LW per unit

%one panel per metric: a faint line per unit (animal x channel x pair, coloured by pair), and the
%median per pair on top. Title: all pairs signrank p and mixed model p
if plot_summary == 1
    colors = lines(npairs);
    figure();
    set(gcf, 'Position', get(0, 'Screensize'));
    %tiledlayout(2,4,'TileSpacing','compact'); %commented out: 10 metrics now
    tiledlayout(2,5,'TileSpacing','compact');
    for k = 1:nmet
        ax = nexttile;
        hold(ax,'on')
        for i = 1:numel(u_pair)
            if isnan(uL(i,k)) == false && isnan(uLW(i,k)) == false
                plot(ax,[1 2],[uL(i,k) uLW(i,k)],'-','Color',[colors(pairnames == u_pair(i),:) .15])
            end
        end
        h = gobjects(npairs,1);
        for j = 1:npairs
            r = (k-1)*npl + j;
            h(j) = plot(ax,[1 2],[median_L(r) median_LW(r)],'-o','Color',colors(j,:),'LineWidth',2.5,'MarkerFaceColor',colors(j,:));
        end
        hold(ax,'off')
        xlim(ax,[.5 2.5]); xticks(ax,[1 2]); xticklabels(ax,[condi_L condi_LW])
        rall = k*npl; %the "all" row
        title(ax,{char(metricnames(k)), sprintf('signrank p = %.3g, lme p = %.3g', p_signrank(rall), lme_p(rall))},'Interpreter','none','FontSize',9)
        if k == 1
            legend(ax,h,pairnames,'Interpreter','none','Location','best')
        end
    end
    sgtitle('LFP autocorrelation, spectrum and correlations: Light vs Light+Whisker, per unit (animal x channel x pair)')
    if plotsave == 1
        if isfolder([save_directory '\Figures_LFPCorr']) == false
            mkdir([save_directory '\Figures_LFPCorr'])
        end
        savefig(gcf,[save_directory '\Figures_LFPCorr\LFPCorr_summary_' regexprep(char(strjoin(pairnames,'_')),'\W','_')],'compact')
    end
end


%% Data Visualization: trial classes (rhythmic response vs P1 response)

%1) per animal: each pair's trial rhythm scores vs its pre-stim windows, and the P1 scores vs the fake no-stim
%   windows, thresholds dashed (check the calls look sensible)
%2) per pair: Venn of rhythmic vs P1-response trials, a row per animal (+ all animals), L and LW side by side
%3) the 4 classes as % of each animal x condition's good trials
if plot_trialclass == 1
    disp('Trial classes: Making some pictures for you')
    if plotsave == 1 && isfolder([save_directory '\Figures_LFPCorr']) == false
        mkdir([save_directory '\Figures_LFPCorr'])
    end
    classcols = [0.85 0.33 0.10; 0.93 0.69 0.13; 0 0.45 0.74; 0.6 0.6 0.6]; %rhythmic+P1, rhythmic only, P1 only, neither

    %1) scores vs their no-stimulus references
    figure();
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(nmice,npairs+1,'TileSpacing','compact','Padding','compact');
    title(tl,'Trial classes: each trial''s score (L black, LW green) vs the no-stimulus reference (grey), threshold dashed','Interpreter','none')
    for m = 1:nmice
        for j = 1:npairs
            ax = nexttile;
            hold(ax,'on')
            if isempty(rhythmnull{m,j}) == false
                rows = trialclass.Animal_Name == mice(m) & trialclass.Pair == pairnames(j);
                ed = linspace(min([rhythmnull{m,j}; trialclass.Rhythm_Score(rows)]), max([rhythmnull{m,j}; trialclass.Rhythm_Score(rows)]), 30);
                histogram(ax,rhythmnull{m,j},ed,'Normalization','probability','FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
                histogram(ax,trialclass.Rhythm_Score(rows & trialclass.LvsLW == condi_L),ed,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0 0 0],'LineWidth',1.5)
                histogram(ax,trialclass.Rhythm_Score(rows & trialclass.LvsLW == condi_LW),ed,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0 0.6 0],'LineWidth',1.5)
                xline(ax,rhythmthr(m,j),'r--','LineWidth',1.5)
            end
            hold(ax,'off')
            title(ax,{char(mice(m)), char("rhythm: " + pairnames(j))},'Interpreter','none','FontSize',8)
            xlabel(ax,'mean z of P2 rhythm vs pre-stim')
        end
        ax = nexttile;
        hold(ax,'on')
        if isempty(p1null{m}) == false
            rows = trialclass.Animal_Name == mice(m);
            ed = linspace(min([p1null{m}; trialclass.P1_Score(rows)]), max([p1null{m}; trialclass.P1_Score(rows)]), 30);
            histogram(ax,p1null{m},ed,'Normalization','probability','FaceColor',[0.7 0.7 0.7],'EdgeColor','none')
            histogram(ax,trialclass.P1_Score(rows & trialclass.LvsLW == condi_L),ed,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0 0 0],'LineWidth',1.5)
            histogram(ax,trialclass.P1_Score(rows & trialclass.LvsLW == condi_LW),ed,'Normalization','probability','DisplayStyle','stairs','EdgeColor',[0 0.6 0],'LineWidth',1.5)
            xline(ax,p1thr(m),'r--','LineWidth',1.5)
        end
        hold(ax,'off')
        title(ax,{char(mice(m)), 'P1 response (all conditions)'},'Interpreter','none','FontSize',8)
        xlabel(ax,'mean P1 template z')
    end
    if plotsave == 1
        savefig(gcf,[save_directory '\Figures_LFPCorr\TrialClass_scores'],'compact')
    end

    %2) Venn: left circle = rhythmic, right circle = P1 response, "neither" underneath
    th = linspace(0,2*pi,100);
    for j = 1:npairs
        venn_animals = [unique(trialclass.Animal_Name(trialclass.Pair == pairnames(j)),'stable'); "all animals"];
        figure();
        set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout(2,numel(venn_animals),'TileSpacing','compact','Padding','compact');
        title(tl,char("Rhythmic response vs P1 response trials: " + pairnames(j)),'Interpreter','none')
        for g = 1:2
            if g == 1
                grp = condi_L;
            else
                grp = condi_LW;
            end
            for a = 1:numel(venn_animals)
                ax = nexttile(tl,(g-1)*numel(venn_animals) + a);
                rows = trialclass.Pair == pairnames(j) & trialclass.LvsLW == grp;
                if venn_animals(a) ~= "all animals"
                    rows = rows & trialclass.Animal_Name == venn_animals(a);
                end
                n = sum(rows);
                nc = zeros(1,4);
                for c = 1:4
                    nc(c) = sum(trialclass.TrialClass(rows) == classnames(c));
                end
                hold(ax,'on')
                patch(ax,-0.5 + cos(th),sin(th),classcols(2,:),'FaceAlpha',0.3,'EdgeColor',classcols(2,:),'LineWidth',1.5) %rhythmic
                patch(ax,0.5 + cos(th),sin(th),classcols(3,:),'FaceAlpha',0.3,'EdgeColor',classcols(3,:),'LineWidth',1.5) %P1 response
                text(ax,-1,0,num2str(nc(2)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold') %rhythmic only
                text(ax,0,0,num2str(nc(1)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold') %both
                text(ax,1,0,num2str(nc(3)),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold') %P1 only
                text(ax,0,-1.35,['neither: ' num2str(nc(4))],'HorizontalAlignment','center','FontSize',10)
                text(ax,-0.9,1.2,'rhythmic','HorizontalAlignment','center','FontSize',9,'Color',classcols(2,:)*0.8)
                text(ax,0.9,1.2,'P1','HorizontalAlignment','center','FontSize',9,'Color',classcols(3,:)*0.8)
                hold(ax,'off')
                axis(ax,'equal')
                axis(ax,'off')
                xlim(ax,[-1.7 1.7])
                ylim(ax,[-1.6 1.4])
                title(ax,{char(venn_animals(a)), char(grp + " (n = " + n + ")")},'Interpreter','none','FontSize',9)
            end
        end
        if plotsave == 1
            savefig(gcf,[save_directory '\Figures_LFPCorr\TrialClass_Venn_' regexprep(char(pairnames(j)),'\W','_')],'compact')
        end
    end

    %3) the 4 classes as % of each animal x condition
    [grpkeys,~,gi] = unique(trialclass(:,["Animal_Name","Condition_Name"]),'stable');
    pct = zeros(height(grpkeys),4);
    for c = 1:4
        pct(:,c) = accumarray(gi, trialclass.TrialClass == classnames(c)) ./ accumarray(gi,1) * 100;
    end
    figure();
    set(gcf, 'Position', get(0, 'Screensize'));
    b = bar(pct,'stacked','EdgeColor','none');
    for c = 1:4
        b(c).FaceColor = classcols(c,:);
    end
    xticks(1:height(grpkeys))
    xticklabels(cellstr(grpkeys.Animal_Name + " " + grpkeys.Condition_Name))
    set(gca,'TickLabelInterpreter','none')
    xtickangle(30)
    ylabel('% of good trials')
    ylim([0 100])
    legend(b,cellstr(classnames),'Location','eastoutside')
    title('Trial classes per animal x condition')
    box off
    if plotsave == 1
        savefig(gcf,[save_directory '\Figures_LFPCorr\TrialClass_bars'],'compact')
    end
end


%% Data Visualization: ranked-shank heatmaps (same as the summary plots)

%the heatmap functions take the per-trial table as is: each cell = mean over trials of the metric.
%They only plot their own condition names (L_4/LW_4..., or L_8/LW_8 with _TTX), so name your conditions that way
if heatplot > 0
    if isempty(rankcols)
        error('No ShankRank_ columns (no shank scores picked/found), the heatmaps need one to align the shanks')
    end
    firstpick = find(rankcols == "ShankRank_SponActivity",1); %spon activity ranking pre-selected
    if isempty(firstpick)
        firstpick = 1;
    end
    [sel,~] = listdlg('ListString',cellstr(rankcols),'SelectionMode','multiple','Name','Shank alignment', ...
        'PromptString','Centre the shanks on the top shank of which label(s)?','ListSize',[300 150], ...
        'InitialValue',firstpick);
    alignlabels = rankcols(sel);

    heatanimals = unique(tableres.Animal_Name,'stable')'; %animals that made it into the table
    heatfigs = cell(1,numel(alignlabels)); %{1 x alignment labels}, one plot function output each
    heatdata = cell(1,numel(alignlabels));
    for k = 1:numel(alignlabels)
        disp(['Heatmaps: centring on ' char(alignlabels(k))])
        if heatplot == 1
            [heatfigs{k}, heatdata{k}] = NNx_plotLWIntensityCompare_multianimal(heatanimals, tableres, region, heat_metrics, ...
                [], true, save_heatpng == 1, alignlabels(k), true, []);
        elseif heatplot == 2
            [heatfigs{k}, heatdata{k}] = NNx_plotLWTTXCompare_multianimal(heatanimals, tableres, region, heat_metrics, ...
                [], true, save_heatpng == 1, alignlabels(k), true, []);
        end
        set(heatfigs{k},'Name',"LFP corr, aligned on " + alignlabels(k)) %figure titles look the same, so label the windows
    end
end


%% Save: the per-trial table as a csv for R

if save_csv == 1
    dlgtitle = 'What file name to use for the data?';
    promt = {'Give the file a name (no spaces please, no .csv)'};
    fieldsize = [1 150];
    definput = {'LFPCorr_LvsLW'};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    filename = inputdlg(promt,dlgtitle,fieldsize,definput,opts);

    %does it exist?
    exists = isfile(fullfile(csv_directory, [cell2mat(filename) '.csv']));
    while exists == 1
        dlgtitle = 'Already made a file called that, try again';
        promt = {'Give the file a name (that doesnt already exist this time)'};
        filename = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
        exists = isfile(fullfile(csv_directory, [cell2mat(filename) '.csv']));
    end

    writetable(tableres,fullfile(csv_directory, [cell2mat(filename) '.csv']))
    writetable(trialclass,fullfile(csv_directory, [cell2mat(filename) '_trialclass.csv'])) %one row per good trial, its rhythm / P1 class
    disp('saved')
end


%% report

disp('===== NNx_LFPCorrelation_LvsLW run report =====')
disp('files used:')
disp(filesused)
disp('good trials per animal (rows) x condition (columns):')
disp(array2table(ntrials_all,'RowNames',cellstr(mice),'VariableNames',condis))
if rhythm_filter == 1
    disp('analyzed (rhythmic) trials per animal (rows) x condition (columns):')
else
    disp('analyzed trials per animal (rows) x condition (columns) (rhythm_filter = 0, all of them):')
end
disp(array2table(ntrials,'RowNames',cellstr(mice),'VariableNames',condis))
disp(['trial classes: rhythmic = P2 rhythm beats ' num2str(rhythm_null_prct) '% of pre-stim windows, P1 response = beats ' ...
    num2str(p1_null_prct) '% of fake no-stim windows'])
disp(stats_class)
if height(stats_class_glme) > 0
    disp('rhythmic or not on the P1-response trials, across animals (logistic mixed model, animal random):')
    disp(stats_class_glme)
end
disp(stats_paired)
disp(stats_animal)
disp(stats_chanchan)
disp(['fraction of trials with a secondary ACF peak: L ' num2str(mean(u_pkfracL,'omitnan'),3) ...
    ', LW ' num2str(mean(u_pkfracLW,'omitnan'),3) ' (mean over units)'])
disp(['fraction of trials whose ACF never drops below 1/e within ' num2str(acf_maxlag_ms) ' ms: ' ...
    num2str(mean(isnan(tableres.P2_ACFDecay)),3)])
%line noise: if it's much bigger in one condition it can leak into the ACF and correlations (it's left out of the spectrum metrics)
linemed = nan(nmice,ncond); %[animals x conditions] median over trials and channels, dB
for m = 1:nmice
    for con = 1:ncond
        if isempty(linedb{m,con}) == false
            linemed(m,con) = median(linedb{m,con},'all');
        end
    end
end
disp(['line noise: median ' num2str(spec_f(iline)) ' Hz spike over its neighbours (dB), animals (rows) x conditions (columns):'])
disp(array2table(round(linemed,1),'RowNames',cellstr(mice),'VariableNames',condis))
for i = 1:numel(run_notes)
    disp(['NOTE: ' char(run_notes(i))])
end
%NOTE: units (and trials within a channel) aren't independent, so treat signrank / ranksum p values as a
%rough check; lme_p is the one that accounts for animal and channel
%NOTE: with rhythm_filter = 1 the trials were picked on their P2 rhythm, so P2_SpecOscPow (the thing they were picked
%on) and the ACF peak height get cut off at the threshold, which pulls L and LW together. The structure numbers
%(P2_SpecPkFreq, P2_ACFPk2Lag, P2_ACFPkDamp, P2_ACFDecay) and the stats_class proportions are the clean read-outs
