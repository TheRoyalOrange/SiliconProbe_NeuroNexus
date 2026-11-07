% ManifoldPCA_LvsLW.m
%
% Description: quick look at the population response in light alone (L) vs light + whisker (LW), for
%   the channels of one probe (region, e.g. V1). Every trial you pick is drawn as a line through PC
%   space (2D or 3D), all conditions in the same space: L trials grey to black and LW trials pink to red
%   as light intensity goes up (the order you list the conditions in). One column per signal. All four
%   signals are made from the LFP files:
%     LFP      - bandpassed LFP
%     CSD      - inverse CSD of that LFP, per trial (CSD_flex.m), more local than the LFP
%     LFP_RMS  - moving RMS envelope of the LFP (how big the signal is, ignoring its sign and phase)
%     CSD_RMS  - moving RMS envelope of the CSD
%   (Spiking is commented out for now, the code is kept for later.)
%   Steps, per animal x signal:
%     1) loads the good trials of each condition, makes the signals, bins them (bin_ms)
%     2) takes the plot_ms stretch of all trials of all conditions, sticks them end to end, z-scores each
%        channel over all of it (so a gain change between conditions isn't normalised away)
%     3) runs one PCA on that (channels = variables, time bins = samples), so every condition shares the
%        same axes. Conditions with more trials count for more in the fit
%     4) top panel: each trial as a thin line, in random drawing order. Bottom panel: same PCs, just each
%        condition's mean, zoomed in, with a dot at stim onset
%     5) measures the path length of each condition's mean over path_ms (added-up distance between time
%        bins, in the pca_dims PCs), at equal trial counts within an animal, and makes a bar plot. Noise
%        makes a path longer, so read these against a pre-stim path (path_ms) of the same length of time.
%        L vs LW at each intensity is tested by shuffling trials between them (equal trials, PCs fixed),
%        shown on the bars as * (q < 0.05, Benjamini-Hochberg within a signal) or ns. L and LW come from
%        separate files, so drift between files also counts as a difference: use a pre-stim path_ms as control.
%        A second bar plot shows all animals together (mean +- SEM, a line per animal), with an all-animal
%        shuffle test, signrank across animals, and a mixed model with animal as a random effect
%     6) CSD RMS means, L vs LW, one panel per light intensity (same PCs, same axes, equal trials per panel)
%     7) run_pr = 1: participation ratio PR = (sum eig)^2 / sum(eig^2), about how many dimensions the
%        variance is spread over (1 = one PC holds it all, nChans = spread evenly), per condition and window
%        (win_names), all channels, equal trials within an animal, from three covariances:
%          Total - all trials end to end. Noise - each trial minus its condition mean.
%          Signal - covariance BETWEEN the means of two random halves of the trials, averaged over splits,
%          so noise averages away instead of counting as signal
%        Then signal vs noise PR compared directly: same-axes plot, ratio plot, signrank and a mixed model across
%        animals (only conditions with SignalFrac >= pr_minfrac), plus an L vs LW shuffle test on signal PR,
%        noise PR and their ratio (pr_test = 1). Signal PR is capped by a window's time bins, noise PR isn't,
%        so be careful comparing windows of different lengths
%
%   OPTIONAL, run_stats = 1: on top of the plots, for each analysis window (win_names), on the
%   trial means (TrialAvg) and on all trials end to end (SingleTrial), this measures in L and LW:
%       participation ratio, as above but no signal/noise split
%       alignment - how much of one condition's variance fits in the other condition's top pc_k
%         PCs, relative to the most it could be (1 = same subspace), both directions averaged
%       principal angle - mean angle between the two top pc_k subspaces (0 deg = same subspace)
%     L and LW always use the same number of trials here, because PR depends on trial count (with unequal
%     counts the real numbers are averaged over n_sub random subsamples). Trials are shuffled between L and
%     LW n_perm times, because two noisy estimates of the SAME manifold still have alignment < 1 and
%     angle > 0, so compare those to their null values, not to 1 and 0. L and LW come from separate
%     files, so drift or brain state between files also counts as a difference: check the baseline window.
%     Also draws trial-mean figures (plot_traj) and paired L vs LW summaries across animals (plot_summary).
%
% Inputs:
%   User-set config (Setup section): mice, region, signal flags, condis_L, condis_LW, plot_ms, pca_dims,
%     bins, toggles
%   Files picked per condition with uipickfiles (all animals at once, matched by file name start),
%   <lfp_dir>\<animal>\<animal>-<stim>-LFP.mat:
%     stim_lfp_stimchunks (double, trials x 15000 x channels) - LFP, uV [inferred: units not stated
%       upstream], 1 kHz, not filtered upstream
%     tr_keep (double, 1 x nTrials) - trial numbers (indices into dim 1)
%     tr_remove (0/1, 1 x nTrials) - mask over tr_keep (updated format from
%       NNx_TrialRemove_Conditional(_callable).m); only tr_keep(tr_remove == 0) are used. Old format
%       (different length, or empty) stops the script with an error
%   <probeinfo_dir>\<animal>-ProbeInfo.mat, ProbeInfo fields:
%     .Areas (cell of char, 1 x probes) - region per probe; region picks the probe
%     .ProbeMaps (cell, 1 x probes) - [rows x shanks] raw channel IDs, 0 = no site
%     .ProbeIds (cell, 1 x probes) - raw channel IDs on the probe, the channel order CSD_flex.m expects.
%       Must be the same channels as ProbeMaps or the script stops with an error (CSD signals only)
%     .Ch_Remove (cell, 1 x probes) - raw channel IDs to leave out (can be shorter than probes, or empty)
%   INFERRED DATA CONTRACTS (not stated in the files, same as LFPCorrelation_LvsLW.m):
%     - stim onset is sample 5001 at 1 kHz (BaseAnalysis pre = 5 s)
%     - channel index (dim 3) of the data = raw channel ID (true when ChanIds = 1:Chans)
%     - Ch_Remove is indexed by probe number
%     - ProbeMaps rows are equally spaced contacts in physical order along the shank (from CSD_flex.m)
%   Not applied: tr_remove_conditional (extra exclusions from NNx_TrialRemove_Conditional_callable.m)
%
% Outputs (workspace; nothing saved unless plotsave / save_csv):
%   signames (string, 1 x signals) - the signals that were run, order of dim 3 of the cells below
%   stimdat (cell, animals x conditions x signals) - (single, trials x nbins x kept channels) good trials,
%     binned, before z-scoring. LFP in uV, CSD in uA/mm^3 (positive = source), RMS in the same units.
%     Files of the same animal + condition are pooled. CSD signals still hold NaN channels here.
%   time_bins (double, 1 x nbins) - bin centres, ms re: stim
%   chanlist (cell, 1 x animals) - (double, kept channels x 1) raw channel IDs, order of dim 3 of stimdat
%   pcscore (cell, animals x signals) - (double, trials x plot bins x pca_dims) PC scores of every trial
%     of every condition that animal has
%   pccond (cell, animals x signals) - (double, trials x 1) condition of each of those trials, index into condis
%   pcexpl (cell, animals x signals) - (double, PCs x 1) percent variance per PC
%   condcol (double, conditions x 3) - RGB colour of each condition, same order as condis
%   pathlen (double, animals x conditions x signals) - path length of each condition's mean over path_ms,
%     z-scored PC units (raw PC units if znorm = 0), NaN if the animal doesn't have that condition. Only
%     compare within an animal x signal: each has its own PCA
%   pathntr (double, animals x 1) - trials per condition behind each animal's path lengths
%   pathres (table, one row per animal x signal x condition) - Animal_Name, DataType, Condition_Name,
%     nTrials_used, PathLength
%   stats_path (table) - L vs LW shuffle test. One row per animal x signal x pair, plus Pair = "All intensities"
%     rows per animal, plus Animal_Name = "All animals" rows per signal x pair and for all intensities (differences
%     summed, these speak only for these animals). Columns: Animal_Name, DataType, Pair (string); nTrials_used
%     (per condition, can be more than the bars use), Path_L, Path_LW (means over animals on "All animals" rows),
%     Path_diff (LW - L, the "All" rows sum these), Null_sd (SD of the difference under shuffled labels), p
%     (two-sided, uncorrected), q (Benjamini-Hochberg within a signal, separately for single-animal and
%     all-animal rows, NaN on "All intensities" rows)
%   stats_path_paired (table, one row per signal x pair) - across animals: nAnimals, medians, animals with LW
%     longer, signrank p (needs 6+ animals to be able to reach p < 0.05)
%   stats_path_lme (table, one row per signal) - fitlme(log(path) ~ Group (+ Pair) + (1 + Group|Animal), REML,
%     Satterthwaite DF) on pathlen: nAnimals, LWoverL (ratio), CI_lo, CI_hi (95%), DF, p. Needs 2+ animals
%   path_lme_models (cell, 1 x signals) - the fitted models
%   Only with run_pr = 1:
%   prmat (double, animals x conditions x signals x windows x 3) - PR (dims) in order prmeas =
%     Total, Signal, Noise. NaN if the animal doesn't have that condition
%   prres (table, one row per animal x signal x window x condition) - Animal_Name, DataType, Window,
%     Condition_Name, Group (L or LW) (string); nChans, nTrials_used, PR_total, PR_signal, PR_noise (dims),
%     SignalFrac (signal variance / total variance, 0-1; PR_signal means little when this is near 0)
%   prfrac (double, animals x conditions x signals x windows) - SignalFrac, as an array
%   stats_pr_signoise (table, one row per signal x window x condition) - across animals, only those with
%     SignalFrac >= pr_minfrac: nAnimals, nEnoughSignal, PR_signal_med, PR_noise_med, nSignalLower, signrank_p
%     (needs 6+ animals to be able to reach p < 0.05)
%   stats_pr_lme (table, one row per signal x window, 2+ animals) - fitlme(log(PR) ~ Measure (+ Group) +
%     (1 + Measure|Animal), REML, Satterthwaite DF), same conditions: nAnimals, nConditionRows, SignalOverNoise
%     (ratio), CI_lo, CI_hi (95%), DF, p
%   pr_lme_models (cell, signals x windows) - the fitted models
%   stats_pr (table, pr_test = 1) - L vs LW shuffle test, rows per animal x signal x window x pair x Measure
%     (Signal, Noise, Signal/Noise), plus Pair = "All intensities" and Animal_Name = "All animals" rows (differences
%     summed). Columns: Animal_Name, DataType, Window, Pair, Measure (string); nTrials_used, Value_L, Value_LW
%     (one-split PR, averaged over picks; means over animals on "All animals" rows), Diff (LW - L), Null_sd,
%     p (two-sided, uncorrected), q (Benjamini-Hochberg within signal x window x measure, separately for
%     single-animal and all-animal rows, NaN on "All intensities" rows)
%   filesused (table) - animal, condition, file, trials in file, trials used
%   run_notes (string) - problems hit along the way, shown in the report at the end
%   figures (plotsave = 1) - <save_directory>\Figures_Manifold\
%   Only with run_stats = 1:
%   stats_win (table, one row per animal x pair x signal x window x space) - columns:
%     Animal_Name, DataType (= signal name), Region, Pair, Window, Space (string)
%     nChans, nTrials_L, nTrials_LW, nTrials_used (double) - nChans = channels used for this signal,
%       nTrials_used = trials per condition
%     PR_L, PR_LW (dims), PR_diff (LW - L), PR_null_sd (SD of PR_diff under shuffled labels),
%       PR_p (two-sided permutation, uncorrected)
%     TotVar_L, TotVar_LW (z-units^2) - total variance, to tell a shape change from a size change
%     Align (0-1), Align_null (median under shuffled labels), Align_p (one-sided: L and LW LESS aligned
%       than chance, uncorrected)
%     Angle (deg), Angle_null, Angle_p (one-sided: bigger than chance, uncorrected)
%   stats_paired (table, one row per signal x pair x window x space) - across animals: n animals, median
%     PR_L and PR_LW, animals with LW lower, signrank p on PR (needs 6+ animals to be able to reach
%     p < 0.05), median Align and Align_null
%   trajL, trajLW (cell, animals x pairs x signals) - (double, plot bins x pc_k) shared-PC scores of the
%     trial means, equal trials in L and LW
%   jointexpl (cell, animals x pairs x signals) - (double, PCs x 1) percent variance per shared PC
%   csv (save_csv = 1) - stats_win, you'll get a file name box
%
% Dependencies: NNx_BaseAnalysis*.m (LFP files, ProbeInfo); NNx_TrialRemove_Conditional(_callable).m
%   (tr_remove mask format); CSD_flex.m (Functional_Connectivity folder, CSD signals); uipickfiles.m;
%   Statistics Toolbox (pca, signrank, fitlme); Signal Processing Toolbox (butter, zp2sos, filtfilt,
%   envelope); MATLAB R2019b+ (tiledlayout).
%

%% Setup: Specify Animals, Conditions, and Analysis Info

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Who is the data? (exactly how their file names start: <animal>-<stim>-...)
mice = {"20260508-p9", "20260509-p10", "20260821-p10"}; % {"20260226-p12", "20260411-p9", "20260423-p12"};

%which probe should I analyze? (must match an entry of ProbeInfo.Areas for every mouse, or a probe number like "1")
region = "V1";

%which signals should I run? (1 = yes, 0 = no)
LFP = 0; %bandpassed LFP. volume conduction makes channels very alike
CSD = 0; %current source density of that LFP, more local
LFPRMS = 1; %RMS envelope of the LFP
CSDRMS = 0; %RMS envelope of the CSD

%spiking is off for now, LFP and CSD only. Kept for later, along with its loading block further down
%datatype = 2;
    % 1 = LFP (volume conduction makes channels very alike, so expect a low PR)
    % 2 = MUA (spikes, more local)

%which conditions are LIGHT (L) and which are LIGHT + WHISKER (LW)? (you'll pick each one's files next)
%same number of each -> paired in list order (1st L with 1st LW, ...) and compared pair by pair.
%otherwise all L vs all LW as one lumped pair (pairing only matters for the optional stats).
%LIST THEM LOWEST INTENSITY FIRST: the order sets the colour, grey to black for L, pink to red for LW
%condis_L = {'L_4','L_8','L_12','L_15'};
%condis_LW = {'LW_4','LW_8','LW_12','LW_15'};
condis_L = {'L_8' 'LW_8'};
condis_LW = {'L_8_TTX' 'LW_8_TXX'};

%where are the data files?
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';
%spike_dir = 'E:\Roy\Processed Silicon Probe Data\Spiking'; %spiking off for now
probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%Where should I, Matlab, put the results?
csv_directory = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\csvfiles_forR';
save_directory = 'E:\Roy\Processed Silicon Probe Data\BundledAnimalData\MultAnimal'; %figures go in \Figures_Manifold

%%%Plotting details%%%
%2D or 3D PCA plot?
pca_dims = 3;
    % 2 = PC1 vs PC2
    % 3 = PC1 vs PC2 vs PC3
plot_ms = [-200 1000]; %ms re: stim, the stretch that gets plotted, and that the PCA is fit to
%colours: each list goes from the first colour at its lowest intensity to the second at its highest
col_grey = [0.7 0.7 0.7];
col_black = [0 0 0];
col_pink = [1 0.7 0.78];
col_red = [0.85 0 0];
%save the figures?
plotsave = 0;

%%%Analysis details%%%
bin_ms = 5; %ms per time bin (whole number). All the PCA runs on these bins

%z-score each channel? (without it, the loudest few channels decide the PCs)
znorm = 1;
    % 0 = no
    % 1 = divide each channel by its SD over L + LW together
%norm_soft = 5; %MUA only, spikes/s added to the SD so near-silent channels don't get blown up. Spiking off for now

%LFP cleanup
lfp_band = [2 100]; %Hz, zero-phase Butterworth bandpass (same as LFPCorrelation_LvsLW.m). A lower highpass rings further back into the baseline
filt_pad_ms = 1000; %ms of extra data on each side for the filter to settle

%RMS envelope
rms_win_ms = 25; %ms, moving RMS window (same as the summary script)

%inverse CSD (CSD_flex.m, same numbers as ITPC_PhaseReset_singleanimal.m)
csd_spacing_um = 100; %distance between adjacent ProbeMaps rows, um
csd_radius_um = 300; %lateral radius of the active column of tissue, um (not the contact size)
csd_conductivity = 0.3; %extracellular conductivity, S/m
csd_min_chans = 3; %shanks with fewer kept channels get no CSD, so their channels are left out of the CSD signals

%MUA cleanup (spiking off for now)
%mua_smooth_ms = 25; %ms, gaussian smoothing window of the binned rate (0 = none)

%%%Path length%%%
path_ms = [0 3000]; %ms re: stim, the stretch of each condition's mean to measure (has to be inside plot_ms)
path_nsub = 100; %random picks of trials to average over, for conditions with more trials than the animal's lowest
path_nperm = 500; %label shuffles for the L vs LW test on the bar plot (lower it if you're in a hurry)

%%%Participation ratio%%%
%want the participation ratio report, split into signal and noise?
run_pr = 1;
    % 0 = no
    % 1 = yes
pr_nsub = 20; %random picks and splits of trials to average over

%which windows should I measure? (also used by the optional stats)
%baseline is your control: L and LW shouldn't differ much before the stim
win_names = ["Baseline","Resp","P1","P2"];
win_ms = [-1500 -500; 75 3000; 75 350; 350 3000]; %ms re: stim, one row per window. Baseline stops well before the stim because the zero-phase LFP filter smears the response backward in time
pr_minfrac = 0.05; %SignalFrac below this means there's too little signal to trust PR_signal: hollow markers, and left out of the signal vs noise tests
%want the L vs LW shuffle test on signal PR, noise PR, and their ratio? (this is the slow part)
pr_test = 1;
    % 0 = no
    % 1 = yes
pr_nperm = 100; %label shuffles for that test

%%%Optional stats%%%
%want alignment, angles, L vs LW shuffles, and their figures too?
run_stats = 1;
    % 0 = no
    % 1 = yes (the settings below only matter if this is 1)

pc_k = 5; %how many top PCs make up the "subspace" for alignment and angles
n_perm = 1000; %label shuffles for the permutation test (lower it if you're in a hurry)
n_sub = 50; %subsamples to average the real L vs LW numbers over when trial counts differ
min_trials = 10; %skip an animal x pair with fewer trials than this in L or LW
rng_seed = 1; %same shuffles every run, so results repeat

%trial-mean trajectories, PC time courses, spectra? (one figure per animal x pair x signal)
plot_traj = 1;
%paired L vs LW across animals? (one PR figure and one alignment figure per pair x signal)
plot_summary = 1;
%save the stats table as a csv for R?
save_csv = 0;
    % 0 = no
    % 1 = yes (you'll get a file name box)

%%%nothing else to fill in%%%


%% Data PreProcessing: line up the conditions into L/LW pairs

rng(rng_seed);
mice = string(mice); %string array, easier to compare
nmice = numel(mice);
condis = [condis_L condis_LW]; %all conditions, L first
ncond = numel(condis);
nwin = numel(win_names);
spacenames = ["TrialAvg","SingleTrial"];

%which signals are on? s below always counts through just these
allsignames = ["LFP","CSD","LFP_RMS","CSD_RMS"];
sigcode = find([LFP CSD LFPRMS CSDRMS] == 1); %1 = LFP, 2 = CSD, 3 = LFP RMS, 4 = CSD RMS
signames = allsignames(sigcode);
nsig = numel(signames);
needCSD = CSD == 1 || CSDRMS == 1; %skip the CSD maths if nobody wants it

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

%spiking off for now, every signal comes from the LFP files
%{
if datatype == 1
    typename = "LFP";
    pick_dir = lfp_dir;
    pick_filter = '-LFP\.mat$';
else
    typename = "MUA";
    pick_dir = spike_dir;
    pick_filter = '-spiking_results\.mat$';
end
%}

run_notes = strings(0,1); %problems for the report at the end

%create the figure folder
if plotsave == 1
    if isfolder([save_directory '\Figures_Manifold']) == false
        mkdir([save_directory '\Figures_Manifold'])
    end
end


%% Data Extraction: pick the files of each condition (all animals at once)

superCondis_dir = cell(1,ncond);
for con = 1:ncond
    superCondis_dir{con} = uipickfiles('FilterSpec',lfp_dir,'Prompt',['Choose ' condis{con} ' LFP files (all animals)'],'REFilter','-LFP\.mat$');
    if iscell(superCondis_dir{con}) == false
        error(['No files picked for ' condis{con}])
    end
end

%which mouse is each file from? (file names start with animal-)
condiFileMice = cell(1,ncond); %per file, index into mice
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
    end
end


%% Data Extraction: each animal's probe and kept channels (ProbeInfo)

pinfos = cell(1,nmice); %each animal's ProbeInfo, CSD_flex needs it
mouseProbe = zeros(1,nmice); %probe index used per animal
chanlist = cell(1,nmice); %column of kept raw channel IDs (sorted)
csdorder = cell(1,nmice); %where each channel in CSD_flex's order sits in chanlist
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
    chremove = [];
    if numel(pinfos{m}.Ch_Remove) >= prb
        chremove = pinfos{m}.Ch_Remove{prb};
    end
    chanlist{m} = setdiff(nonzeros(pinfos{m}.ProbeMaps{prb}), chremove(:)); %sorted column of raw channel IDs

    %CSD_flex wants channels in ProbeIds order, so work out how to get there and back
    if needCSD == true
        pids = pinfos{m}.ProbeIds{prb}(:);
        pids = pids(ismember(pids,chremove(:)) == false); %ProbeIds minus Ch_Remove, CSD_flex's order
        if isequal(sort(pids),chanlist{m}) == false
            error([char(mice(m)) ': ProbeIds and ProbeMaps of probe ' num2str(prb) ' list different channels, so I can''t run CSD_flex'])
        end
        [~,csdorder{m}] = ismember(pids,chanlist{m});
    end
end
clear pinfo


%% Data Extraction: load the good trials of every file, make the signals, and bin them

%time bins kept (the plot stretch, plus the PR/stats windows if you want those), plus padding so filter/envelope edges stay outside
keep_ms = plot_ms;
if run_pr == 1 || run_stats == 1
    keep_ms = [min([win_ms(:); plot_ms(:)]) max([win_ms(:); plot_ms(:)])];
end
nbins = ceil((keep_ms(2)-keep_ms(1))/bin_ms); %bins kept
time_bins = keep_ms(1) + bin_ms*((1:nbins)-0.5); %bin centres, ms re: stim

padbins = ceil(filt_pad_ms/bin_ms);
[z,p,k] = butter(2, lfp_band/(1000/2), 'bandpass'); %2nd order per edge, as the summary script
[lfp_sos,lfp_g] = zp2sos(z,p,k); %second-order sections, stable at low cutoffs
%spiking off for now
%padbins = ceil(2*mua_smooth_ms/bin_ms);
nbins_pad = nbins + 2*padbins; %bins including padding
start_ms = keep_ms(1) - padbins*bin_ms; %ms re: stim of the first loaded sample
inkeep = padbins + (1:nbins); %where the kept bins sit inside the padded ones
perbin = bin_ms; %samples per bin, 1 kHz
samps = 5001 + start_ms + (0:(nbins_pad*perbin-1)); %5001 = stim onset, change if the epochs change

stimdat = cell(nmice,ncond,nsig); %{animal x condition x signal} [trials x nbins x kept channels], single
filesused = table();

tic
for con = 1:ncond
    for file = 1:numel(superCondis_dir{con})
        m = condiFileMice{con}(file);
        fpath = superCondis_dir{con}{file};
        [~,fname] = fileparts(fpath);
        disp(['LFP: Loading ' condis{con} ', ' fname])
        chs = chanlist{m};

        %which trials are good? tr_remove has to be the 0/1 mask over tr_keep
        d = load(fpath,'tr_keep','tr_remove');
        if isfield(d,'tr_remove') == false || numel(d.tr_remove) ~= numel(d.tr_keep)
            error([fname ': tr_remove is not a 0/1 mask over tr_keep (old format). Run the file through NNx_TrialRemove_Conditional first'])
        end
        tr_keep = d.tr_keep(:)';
        trs = tr_keep(d.tr_remove(:)' == 0); %good trials

        %LFP files aren't v7.3, so load it all and take what we need
        d = load(fpath,'stim_lfp_stimchunks');
        dat = d.stim_lfp_stimchunks(trs,samps,chs); %[good trials x samples x kept channels]
        clear d
        [ntr,nsamp,nch] = size(dat);
        X = reshape(permute(dat,[2 1 3]), nsamp, []); %[samples x (trials*channels)], filtfilt works down columns
        X = filtfilt(lfp_sos, lfp_g, X); %zero phase
        lfpflt = permute(reshape(X, nsamp, ntr, nch), [2 1 3]); %back to [trials x samples x channels], uV
        clear dat

        %inverse CSD of the bandpassed LFP, trial by trial and shank by shank. NaN on shanks with too few channels
        if needCSD == true
            csdflt = nan(ntr,nsamp,nch);
            csdflt(:,:,csdorder{m}) = CSD_flex(lfpflt(:,:,csdorder{m}), pinfos{m}, mouseProbe(m), 'ChRemove', true, ...
                'spacing', csd_spacing_um, 'radius', csd_radius_um, 'conductivity', csd_conductivity, 'minchans', csd_min_chans);
        end

        %make each signal that's switched on, then bin it and pool it
        for s = 1:nsig
            switch sigcode(s)
                case 1 %LFP
                    sigdat = lfpflt;
                case 2 %CSD
                    sigdat = csdflt;
                case 3 %LFP RMS
                    Y = envelope(X, rms_win_ms, 'rms'); %X is still the bandpassed LFP as [samples x (trials*channels)]
                    sigdat = permute(reshape(Y, nsamp, ntr, nch), [2 1 3]);
                case 4 %CSD RMS
                    csdcols = reshape(permute(csdflt,[2 1 3]), nsamp, []); %[samples x (trials*channels)]
                    hascsd = all(isfinite(csdcols),1); %envelope can't take NaN, so skip the channels with no CSD
                    Y = nan(size(csdcols));
                    Y(:,hascsd) = envelope(csdcols(:,hascsd), rms_win_ms, 'rms');
                    sigdat = permute(reshape(Y, nsamp, ntr, nch), [2 1 3]);
                    clear csdcols
            end
            binned = reshape(mean(reshape(single(sigdat),ntr,perbin,nbins_pad,nch),2),ntr,nbins_pad,nch); %mean per bin
            stimdat{m,con,s} = cat(1, stimdat{m,con,s}, binned(:,inkeep,:)); %pool the files of the same animal + condition
        end
        clear X Y lfpflt csdflt sigdat binned

        %spiking off for now, every signal above comes from the LFP. Kept for later
        %{
        %spike files are huge (450000 samples), so read only our stretch.
        %matfile can only read evenly spaced ranges: read the block, then pick from it in memory
        perbin = bin_ms*30; %samples per bin, 30 kHz
        samp1 = 150001 + start_ms*30; %150001 = stim onset, change if the epochs change
        tr_rng = min(trs):max(trs);
        ch_rng = min(chs):max(chs);
        sdat = matfile(fpath);
        dat = sdat.stim_spike_stimchunks(tr_rng,samp1:(samp1+nbins_pad*perbin-1),ch_rng);
        dat = dat(trs-tr_rng(1)+1,:,chs-ch_rng(1)+1); %[good trials x samples x kept channels]
        [ntr,~,nch] = size(dat);
        binned = reshape(sum(reshape(single(dat),ntr,perbin,nbins_pad,nch),2),ntr,nbins_pad,nch)/(bin_ms/1000); %spikes/s per bin
        if mua_smooth_ms > 0
            binned = smoothdata(binned,2,'gaussian',mua_smooth_ms/bin_ms); %smooth along time within each trial
        end
        %}

        filesused = [filesused; table(mice(m),string(condis{con}),string(fname),numel(tr_keep),numel(trs), ...
            'VariableNames',{'Animal_Name','Condition_Name','File','nTrials_kept','nTrials_used'})]; %#ok<AGROW>
    end
end
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data Visualization: PCA of all trials, coloured by L or LW

%one figure per animal, one column per signal. Every condition you picked goes into one PCA, so they all
%share the same axes. Each trial is one line through PC space
inplot = find(time_bins >= plot_ms(1) & time_bins < plot_ms(2)); %bins to show
tplot = time_bins(inplot); %ms
[~,onset_bin] = min(abs(tplot)); %bin closest to stim onset, marked on the means
col_L = [0 0 0]; %black = light (optional stats figures)
col_LW = [0.85 0.15 0.15]; %red = light + whisker (optional stats figures)
trial_alpha = 0.3; %how see-through the single trial lines are (0 = invisible, 1 = solid)

%colour of each condition: L goes grey to black, LW goes pink to red, in the order you listed them.
%just one condition in a list gets the darkest colour
condcol = zeros(ncond,3); %[conditions x RGB], same order as condis
fracL = linspace(0,1,numel(condis_L)); %0 = lowest intensity colour, 1 = highest
for j = 1:numel(condis_L)
    condcol(j,:) = col_grey + fracL(j)*(col_black - col_grey);
end
fracLW = linspace(0,1,numel(condis_LW));
for j = 1:numel(condis_LW)
    condcol(numel(condis_L)+j,:) = col_pink + fracLW(j)*(col_red - col_pink);
end

pcscore = cell(nmice,nsig); %{animal x signal} [trials x plot bins x pca_dims]
pccond = cell(nmice,nsig); %{animal x signal} [trials x 1] which condition each trial is, index into condis
pcexpl = cell(nmice,nsig); %{animal x signal} percent variance per PC

disp('PCA: Making some pictures for you')
for m = 1:nmice
    %which conditions does this animal have?
    havecond = [];
    for con = 1:ncond
        if isempty(stimdat{m,con,1}) == false
            havecond(end+1) = con; %#ok<SAGROW>
        end
    end
    if isempty(havecond) == true
        run_notes(end+1) = mice(m) + ": no PCA plot, no files picked"; %#ok<SAGROW>
        continue
    end

    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tiledlayout(2,nsig) %top row = all trials, bottom row = means, one column per signal

    for s = 1:nsig
        %stick all of this animal's conditions together, just the plot stretch
        pooled = [];
        trcond = [];
        for con = havecond
            pooled = cat(1, pooled, stimdat{m,con,s}(:,inplot,:)); %[all trials x plot bins x channels]
            trcond = [trcond; repmat(con,size(stimdat{m,con,s},1),1)]; %#ok<AGROW>
        end

        %channels on shanks with too few kept channels have no CSD, so they sit this signal out
        chgood = all(isfinite(pooled),[1 2]);
        if all(chgood) == false
            run_notes(end+1) = mice(m) + ", " + signames(s) + ": left out " + nnz(chgood == false) + " channels with no data (too few kept channels on their shank)"; %#ok<SAGROW>
            pooled = pooled(:,:,chgood);
        end
        [ntr,nb,nch] = size(pooled);

        %all trials end to end, then one PCA, so every condition shares the same axes
        Xall = double(reshape(pooled,[],nch)); %[(trials*bins) x channels]
        if znorm == 1
            Xall = Xall ./ std(Xall,0,1); %z-score over all conditions together. PCA takes out the mean
        end
        [~,score,latent] = pca(Xall,'NumComponents',pca_dims); %score is [(trials*bins) x pca_dims]
        pcscore{m,s} = reshape(score,ntr,nb,pca_dims); %[trials x plot bins x pca_dims]
        pccond{m,s} = trcond;
        pcexpl{m,s} = 100*latent/sum(latent); %percent variance per PC

        %top panel: every trial, thin and see-through. Random drawing order so no colour always sits on top
        nexttile(s)
        hold on
        for tr = randperm(ntr)
            trcol = [condcol(trcond(tr),:) trial_alpha];
            trialpc = reshape(pcscore{m,s}(tr,:,:),nb,pca_dims); %[plot bins x pca_dims]
            if pca_dims == 3
                plot3(trialpc(:,1),trialpc(:,2),trialpc(:,3),'Color',trcol,'LineWidth',0.5)
            else
                plot(trialpc(:,1),trialpc(:,2),'Color',trcol,'LineWidth',0.5)
            end
        end

        %invisible lines, just so the legend has solid colours
        hcon = gobjects(1,numel(havecond));
        labs = cell(1,numel(havecond));
        for j = 1:numel(havecond)
            hcon(j) = plot(nan,nan,'Color',condcol(havecond(j),:),'LineWidth',2);
            labs{j} = [condis{havecond(j)} ', ' num2str(sum(trcond == havecond(j))) ' trials'];
        end
        hold off
        legend(hcon,labs,'Location','best','Interpreter','none')
        title([char(signames(s)) ', ' num2str(nch) ' channels, all trials'],'Interpreter','none')

        %bottom panel: same PCs, just each condition's mean, zoomed in (it's tiny next to the single trials)
        nexttile(nsig+s)
        hold on
        for con = havecond
            meancon = reshape(mean(pcscore{m,s}(trcond == con,:,:),1),nb,pca_dims); %[plot bins x pca_dims]
            if pca_dims == 3
                plot3(meancon(:,1),meancon(:,2),meancon(:,3),'Color',condcol(con,:),'LineWidth',1.5)
                plot3(meancon(onset_bin,1),meancon(onset_bin,2),meancon(onset_bin,3),'o','Color',condcol(con,:),'MarkerFaceColor',condcol(con,:))
            else
                plot(meancon(:,1),meancon(:,2),'Color',condcol(con,:),'LineWidth',1.5)
                plot(meancon(onset_bin,1),meancon(onset_bin,2),'o','Color',condcol(con,:),'MarkerFaceColor',condcol(con,:))
            end
        end
        hold off
        title([char(signames(s)) ', mean of each condition'],'Interpreter','none')

        %same axis labels on both panels
        for panel = [s nsig+s]
            nexttile(panel)
            xlabel(['PC1 (' num2str(pcexpl{m,s}(1),'%.1f') '%)'])
            ylabel(['PC2 (' num2str(pcexpl{m,s}(2),'%.1f') '%)'])
            if pca_dims == 3
                zlabel(['PC3 (' num2str(pcexpl{m,s}(3),'%.1f') '%)'])
                view(3)
                grid on
            end
        end
    end
    sgtitle([char(mice(m)) ' ' char(region) ': L grey to black, LW pink to red, low to high intensity, dot = stim onset'],'Interpreter','none')

    if plotsave == 1
        fname = sprintf([char(mice(m)) '-' 'trialPCA']);
        savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
    end
end
clear pooled Xall score


%% Data Analysis: path length of each condition's mean trajectory

%path length = how far the condition's mean travels through PC space: add up the straight-line distance
%from each time bin to the next, using the pca_dims PCs you see in the plots.
%a mean made from fewer trials is noisier, wigglier, and so longer, whatever the brain did. So within each
%animal, every condition uses the same number of trials: its fewest-trial condition's count. Conditions
%with more trials than that are averaged over path_nsub random picks of trials
inpath = find(tplot >= path_ms(1) & tplot < path_ms(2)); %which of the plotted bins to measure
pathlen = nan(nmice,ncond,nsig); %[animals x conditions x signals] path length, z-scored PC units
pathntr = nan(nmice,1); %trials per condition behind each animal's path lengths

for m = 1:nmice
    for s = 1:nsig
        if isempty(pcscore{m,s}) == true
            continue %no files for this animal
        end
        trcond = pccond{m,s};
        havecond = unique(trcond)'; %conditions this animal has

        %fewest trials in any of this animal's conditions
        nuse = numel(trcond);
        for con = havecond
            nuse = min(nuse,sum(trcond == con));
        end
        pathntr(m) = nuse;

        for con = havecond
            idxcon = find(trcond == con); %this condition's trials
            nsub_use = path_nsub;
            if numel(idxcon) == nuse
                nsub_use = 1; %already the right count, no need to repeat
            end
            lens = zeros(nsub_use,1);
            for it = 1:nsub_use
                pick = idxcon(randperm(numel(idxcon),nuse)); %nuse of this condition's trials
                meancon = reshape(mean(pcscore{m,s}(pick,inpath,:),1),numel(inpath),pca_dims); %[path bins x pca_dims]
                lens(it) = sum(sqrt(sum(diff(meancon,1,1).^2,2))); %distance from each bin to the next, added up
            end
            pathlen(m,con,s) = mean(lens);
        end
    end
end

%same numbers as a table, one row per animal x signal x condition
pathres = table();
for m = 1:nmice
    for s = 1:nsig
        for con = 1:ncond
            if isnan(pathlen(m,con,s)) == false
                pathres = [pathres; table(mice(m),signames(s),string(condis{con}),pathntr(m),pathlen(m,con,s), ...
                    'VariableNames',{'Animal_Name','DataType','Condition_Name','nTrials_used','PathLength'})]; %#ok<AGROW>
            end
        end
    end
end


%% Data Analysis: is the path length different in L vs LW? (shuffle test)

%for each L/LW pair: is LW's path longer or shorter than L's, by more than you'd get from which trials
%happened to land in which condition? Shuffle the pair's trials between L and LW path_nperm times and see how
%often the shuffled difference is at least as big as the real one. The PCs stay the same, shuffling only
%changes which trials go into each mean. L and LW use the same number of trials, so noise is the same on
%both sides. This uses all of the pair's trials, so it can use more than the bars do when the animal's
%lowest condition is lower than this pair's.
%L and LW come from separate files, so drift or brain state between files also looks like a difference. Run
%it again with path_ms before the stim (e.g. [-200 0]) as a control: that should NOT come out significant
nbpath = numel(inpath); %bins measured
stats_path = table();
nullall = zeros(path_nperm,nsig); %[shuffles x signals] LW - L summed over everything, for the "All" rows
realall = zeros(1,nsig);
nullpair = zeros(path_nperm,npairs,nsig); %[shuffles x pairs x signals] LW - L summed over animals, per intensity
realpair = zeros(npairs,nsig);
lenpair = zeros(npairs,2,nsig); %[pairs x L/LW x signals] path lengths summed over animals
npairanimals = zeros(npairs,nsig); %how many animals have each pair
tic
for s = 1:nsig
    for m = 1:nmice
        if isempty(pcscore{m,s}) == true
            continue %no files for this animal
        end
        disp([char(signames(s)) ': Path length shuffle test for ' char(mice(m))])
        trcond = pccond{m,s};
        scores = pcscore{m,s}(:,inpath,:); %[trials x path bins x pca_dims]
        nullanimal = zeros(path_nperm,1); %LW - L summed over this animal's intensities
        realanimal = 0;

        for pr = 1:npairs
            idxL = find(ismember(trcond,grpconds{pr,1})); %this pair's L trials
            idxLW = find(ismember(trcond,grpconds{pr,2})); %this pair's LW trials
            if isempty(idxL) == true || isempty(idxLW) == true
                continue %this animal doesn't have both
            end
            nuse = min(numel(idxL),numel(idxLW)); %both get this many trials
            idxboth = [idxL; idxLW];

            %real L and LW, over random picks of equal trials (one go is enough if the counts already match)
            nsub_use = path_nsub;
            if numel(idxL) == numel(idxLW)
                nsub_use = 1;
            end
            lensreal = zeros(nsub_use,2); %[picks x L/LW]
            for it = 1:nsub_use
                pickL = idxL(randperm(numel(idxL),nuse));
                pickLW = idxLW(randperm(numel(idxLW),nuse));
                meanA = reshape(mean(scores(pickL,:,:),1),nbpath,pca_dims); %[path bins x pca_dims]
                meanB = reshape(mean(scores(pickLW,:,:),1),nbpath,pca_dims);
                lensreal(it,1) = sum(sqrt(sum(diff(meanA,1,1).^2,2)));
                lensreal(it,2) = sum(sqrt(sum(diff(meanB,1,1).^2,2)));
            end
            realdiff = mean(lensreal(:,2) - lensreal(:,1)); %LW - L

            %same thing with shuffled labels
            nulldiff = zeros(path_nperm,1);
            for it = 1:path_nperm
                shuf = idxboth(randperm(numel(idxboth)));
                meanA = reshape(mean(scores(shuf(1:nuse),:,:),1),nbpath,pca_dims);
                meanB = reshape(mean(scores(shuf((nuse+1):(2*nuse)),:,:),1),nbpath,pca_dims);
                nulldiff(it) = sum(sqrt(sum(diff(meanB,1,1).^2,2))) - sum(sqrt(sum(diff(meanA,1,1).^2,2)));
            end

            stats_path = [stats_path; table(mice(m),signames(s),pairnames(pr),nuse,mean(lensreal(:,1)),mean(lensreal(:,2)),realdiff,std(nulldiff), ...
                (sum(abs(nulldiff) >= abs(realdiff))+1)/(path_nperm+1),NaN, ...
                'VariableNames',{'Animal_Name','DataType','Pair','nTrials_used','Path_L','Path_LW','Path_diff','Null_sd','p','q'})]; %#ok<AGROW>
            nullanimal = nullanimal + nulldiff;
            realanimal = realanimal + realdiff;
            nullpair(:,pr,s) = nullpair(:,pr,s) + nulldiff;
            realpair(pr,s) = realpair(pr,s) + realdiff;
            lenpair(pr,:,s) = lenpair(pr,:,s) + mean(lensreal,1);
            npairanimals(pr,s) = npairanimals(pr,s) + 1;
        end

        %this animal, all intensities together: does LW run longer or shorter overall?
        stats_path = [stats_path; table(mice(m),signames(s),"All intensities",NaN,NaN,NaN,realanimal,std(nullanimal), ...
            (sum(abs(nullanimal) >= abs(realanimal))+1)/(path_nperm+1),NaN, ...
            'VariableNames',{'Animal_Name','DataType','Pair','nTrials_used','Path_L','Path_LW','Path_diff','Null_sd','p','q'})]; %#ok<AGROW>
        nullall(:,s) = nullall(:,s) + nullanimal;
        realall(s) = realall(s) + realanimal;
    end

    %every animal together, one intensity at a time. Path_L and Path_LW are the means over animals, Path_diff is
    %the sum. This only speaks for these animals, since it shuffles trials and not animals
    for pr = 1:npairs
        if npairanimals(pr,s) == 0
            continue %nobody has this pair
        end
        stats_path = [stats_path; table("All animals",signames(s),pairnames(pr),NaN,lenpair(pr,1,s)/npairanimals(pr,s),lenpair(pr,2,s)/npairanimals(pr,s),realpair(pr,s),std(nullpair(:,pr,s)), ...
            (sum(abs(nullpair(:,pr,s)) >= abs(realpair(pr,s)))+1)/(path_nperm+1),NaN, ...
            'VariableNames',{'Animal_Name','DataType','Pair','nTrials_used','Path_L','Path_LW','Path_diff','Null_sd','p','q'})]; %#ok<AGROW>
    end

    %every animal together, every intensity together
    stats_path = [stats_path; table("All animals",signames(s),"All intensities",NaN,NaN,NaN,realall(s),std(nullall(:,s)), ...
        (sum(abs(nullall(:,s)) >= abs(realall(s)))+1)/(path_nperm+1),NaN, ...
        'VariableNames',{'Animal_Name','DataType','Pair','nTrials_used','Path_L','Path_LW','Path_diff','Null_sd','p','q'})]; %#ok<AGROW>
end

%correct for doing lots of tests: Benjamini-Hochberg q within each signal, done separately for the single
%animal rows and the all-animal rows. The stars on the charts use q, and the "All intensities" rows keep their plain p
for s = 1:nsig
    for scope = 1:2
        if scope == 1
            rowsq = find(stats_path.DataType == signames(s) & stats_path.Pair ~= "All intensities" & stats_path.Animal_Name ~= "All animals");
        else
            rowsq = find(stats_path.DataType == signames(s) & stats_path.Pair ~= "All intensities" & stats_path.Animal_Name == "All animals");
        end
        [psort,ord] = sort(stats_path.p(rowsq));
        nq = numel(psort);
        qsort = psort.*nq./(1:nq)';
        for i = (nq-1):-1:1
            qsort(i) = min(qsort(i),qsort(i+1)); %q can't go down as p goes up
        end
        stats_path.q(rowsq(ord)) = min(qsort,1);
    end
end

%across animals, one number per animal x pair. Signrank needs 6+ animals to be able to reach p < 0.05
stats_path_paired = table();
for s = 1:nsig
    for pr = 1:npairs
        r = stats_path(stats_path.DataType == signames(s) & stats_path.Pair == pairnames(pr) & stats_path.Animal_Name ~= "All animals",:);
        if height(r) == 0
            continue
        end
        stats_path_paired = [stats_path_paired; table(signames(s),pairnames(pr),height(r),median(r.Path_L),median(r.Path_LW),sum(r.Path_LW > r.Path_L),signrank(r.Path_L,r.Path_LW), ...
            'VariableNames',{'DataType','Pair','nAnimals','Path_L_med','Path_LW_med','nAnimals_LWlonger','signrank_p'})]; %#ok<AGROW>
    end
end
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data Visualization: bar plot of path lengths

%one row per animal, one column per signal, same colours as the PCA plots, L and LW side by side at each
%intensity. Above each pair: * = LW differs from L in the shuffle test (q < 0.05), ns = it doesn't.
%each animal and signal has its own PCA, so only compare bars within a panel

%where each bar goes: pairs side by side, with a gap between pairs
barorder = []; %which condition each bar is
barx = []; %where it sits on the x axis
pos = 0;
for pr = 1:npairs
    conspair = [grpconds{pr,1} grpconds{pr,2}];
    barorder = [barorder conspair]; %#ok<AGROW>
    barx = [barx pos+(1:numel(conspair))]; %#ok<AGROW>
    pos = pos + numel(conspair) + 0.75;
end

figure()
set(gcf, 'Position', get(0, 'Screensize'));
tiledlayout(nmice,nsig)
for m = 1:nmice
    for s = 1:nsig
        nexttile
        hb = bar(barx,pathlen(m,barorder,s),'FaceColor','flat','EdgeColor','none');
        hb.CData = condcol(barorder,:);
        xticks(barx)
        xticklabels(condis(barorder))
        set(gca,'TickLabelInterpreter','none')
        xtickangle(45)
        ylabel('Path length (PC units)')
        title([char(mice(m)) ', ' char(signames(s)) ', ' num2str(pathntr(m)) ' trials each'],'Interpreter','none')

        %bracket and result over each pair
        barmax = max(pathlen(m,barorder,s));
        hold on
        for pr = 1:npairs
            r = stats_path(stats_path.Animal_Name == mice(m) & stats_path.DataType == signames(s) & stats_path.Pair == pairnames(pr),:);
            if height(r) == 0
                continue %this animal doesn't have both
            end
            xpair = barx(ismember(barorder,[grpconds{pr,1} grpconds{pr,2}])); %this pair's bars
            ypair = max(pathlen(m,[grpconds{pr,1} grpconds{pr,2}],s)) + 0.04*barmax; %just above them
            plot([xpair(1) xpair(end)],[ypair ypair],'k')
            if r.q < 0.05
                sigtxt = '*';
            else
                sigtxt = 'ns';
            end
            text(mean(xpair),ypair,sigtxt,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9)
        end
        hold off
        ylim([0 1.25*barmax])
    end
end
sgtitle([char(region) ': path length of each condition mean, ' num2str(path_ms(1)) ' to ' num2str(path_ms(2)) ' ms, first ' num2str(pca_dims) ' PCs. L vs LW shuffle test: * = q < 0.05'])

if plotsave == 1
    fname = sprintf([char(region) '-' 'pathlength']);
    savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
end


%% Data Analysis: mixed model of path length across animals

%does LW change path length overall, treating animals as a sample of mice? Each animal gets its own baseline and
%its own LW effect, as in LFPCorrelation_LvsLW.m: log(path) ~ Group (+ Pair) + (1 + Group|Animal). Log, so the
%effect is a ratio (LW/L), and animals with bigger PC units don't count for more. Uses the bar plot's numbers.
%with few animals this is shaky: there's not much to estimate each animal's own LW effect from
stats_path_lme = table();
path_lme_models = cell(1,nsig);
for s = 1:nsig
    %one row per animal x condition
    lmeAnimal = strings(0,1);
    lmeGroup = strings(0,1);
    lmePair = strings(0,1);
    lmePath = zeros(0,1);
    for m = 1:nmice
        for pr = 1:npairs
            for g = 1:2
                for con = grpconds{pr,g}
                    if isnan(pathlen(m,con,s)) == false
                        lmeAnimal(end+1,1) = mice(m); %#ok<SAGROW>
                        lmeGroup(end+1,1) = string(g); %#ok<SAGROW> 1 = L, 2 = LW
                        lmePair(end+1,1) = pairnames(pr); %#ok<SAGROW>
                        lmePath(end+1,1) = log(pathlen(m,con,s)); %#ok<SAGROW>
                    end
                end
            end
        end
    end

    %needs at least 2 animals with both L and LW to mean anything
    nanimals_lme = numel(unique(lmeAnimal));
    if nanimals_lme < 2
        run_notes(end+1) = signames(s) + ": no path length mixed model, needs at least 2 animals"; %#ok<SAGROW>
        continue
    end
    tbl = table(categorical(lmeAnimal),categorical(lmeGroup),categorical(lmePair),lmePath,'VariableNames',{'Animal','Group','Pair','logPath'});
    if npairs > 1
        path_lme_models{s} = fitlme(tbl,'logPath ~ Group + Pair + (1 + Group|Animal)');
    else
        path_lme_models{s} = fitlme(tbl,'logPath ~ Group + (1 + Group|Animal)');
    end
    [~,~,fx] = fixedEffects(path_lme_models{s},'DFMethod','satterthwaite');
    rowLW = strcmp(fx.Name,'Group_2'); %the LW effect
    stats_path_lme = [stats_path_lme; table(signames(s),nanimals_lme,exp(fx.Estimate(rowLW)),exp(fx.Lower(rowLW)),exp(fx.Upper(rowLW)),fx.DF(rowLW),fx.pValue(rowLW), ...
        'VariableNames',{'DataType','nAnimals','LWoverL','CI_lo','CI_hi','DF','p'})]; %#ok<AGROW>
end


%% Data Visualization: mean path length of all animals together

%one panel per signal. Bars = mean over animals, error bars = SEM, grey lines = each animal's L to LW.
%each animal has its own PCA and channels, so the units aren't exactly the same between animals: the lines
%matter more than the bar heights. Bracket text: q = shuffle test over all animals (trials shuffled, so it speaks
%only for these animals). sr = signrank across animals (animals are the n, needs 6+ to reach p < 0.05).
%panel title: the same two kinds of test, for all intensities together, with the mixed model as the across-animal one
if numel(unique(pathntr(isnan(pathntr) == false))) > 1
    run_notes(end+1) = "Path length: animals have different trial counts (" + strjoin(string(pathntr'),", ") + "). Fewer trials means longer paths, so bar heights between animals aren't comparable, but each animal's L vs LW still is";
end

figure()
set(gcf, 'Position', get(0, 'Screensize'));
tiledlayout('flow')
for s = 1:nsig
    nexttile
    pathvals = pathlen(:,barorder,s); %[animals x bars]
    pathmean = mean(pathvals,1,'omitnan');
    pathsem = std(pathvals,0,1,'omitnan')./sqrt(sum(isnan(pathvals) == false,1));

    hb = bar(barx,pathmean,'FaceColor','flat','EdgeColor','none');
    hb.CData = condcol(barorder,:);
    hold on
    errorbar(barx,pathmean,pathsem,'k','LineStyle','none','LineWidth',1)

    %each animal's L to LW, and the stats above each pair
    barmax = max(pathvals(:));
    for pr = 1:npairs
        barsp = find(ismember(barorder,[grpconds{pr,1} grpconds{pr,2}])); %this pair's bars
        for m = 1:nmice
            plot(barx(barsp),pathvals(m,barsp),'-o','Color',[0.5 0.5 0.5],'MarkerSize',3,'MarkerFaceColor',[0.5 0.5 0.5])
        end

        rq = stats_path(stats_path.Animal_Name == "All animals" & stats_path.DataType == signames(s) & stats_path.Pair == pairnames(pr),:);
        rs = stats_path_paired(stats_path_paired.DataType == signames(s) & stats_path_paired.Pair == pairnames(pr),:);
        if height(rq) == 0
            continue %nobody has this pair
        end
        ypair = max(pathvals(:,barsp),[],'all') + 0.04*barmax; %just above them
        plot([barx(barsp(1)) barx(barsp(end))],[ypair ypair],'k')
        sigtxt = ['q = ' num2str(rq.q,'%.3f')];
        if rq.q < 0.05
            sigtxt = ['* ' sigtxt];
        end
        text(mean(barx(barsp)),ypair,{sigtxt,['sr p = ' num2str(rs.signrank_p,'%.2f')]},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    end
    hold off
    xticks(barx)
    xticklabels(condis(barorder))
    set(gca,'TickLabelInterpreter','none')
    xtickangle(45)
    ylabel('Path length (PC units)')
    ylim([0 1.3*barmax])

    %all intensities together, in the title
    rall = stats_path(stats_path.Animal_Name == "All animals" & stats_path.DataType == signames(s) & stats_path.Pair == "All intensities",:);
    titletxt = {[char(signames(s)) ', all intensities: shuffle p = ' num2str(rall.p,'%.3f')]};
    if height(stats_path_lme) > 0 && any(stats_path_lme.DataType == signames(s))
        rl = stats_path_lme(stats_path_lme.DataType == signames(s),:);
        titletxt{2} = ['mixed model LW/L = ' num2str(rl.LWoverL,'%.2f') ' [' num2str(rl.CI_lo,'%.2f') ' ' num2str(rl.CI_hi,'%.2f') '], p = ' num2str(rl.p,'%.3f')];
    end
    title(titletxt,'Interpreter','none')
end
sgtitle([char(region) ': mean path length over ' num2str(nmice) ' animals, ' num2str(path_ms(1)) ' to ' num2str(path_ms(2)) ' ms, first ' num2str(pca_dims) ' PCs, mean +- SEM, grey = animals'])

if plotsave == 1
    fname = sprintf([char(region) '-' 'pathlength_allanimals']);
    savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
end


%% Data Visualization: CSD RMS means, L vs LW, one panel per light intensity

%same PCs as the PCA plots, so the panels can be compared. One figure per animal, one panel per L/LW pair.
%L and LW in a panel use the same number of trials (a random pick), so neither looks smoother from having more
csdrms_s = find(signames == "CSD_RMS"); %which signal is CSD RMS
if isempty(csdrms_s) == true
    run_notes(end+1) = "No CSD RMS by intensity figure, CSDRMS is switched off";
else
    for m = 1:nmice
        if isempty(pcscore{m,csdrms_s}) == true
            continue %no files for this animal
        end
        trcond = pccond{m,csdrms_s};
        nb = numel(inplot);

        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        tiledlayout('flow') %four intensities makes this 2 x 2
        lims = [inf -inf]; %biggest range of any mean, so every panel gets the same axes
        pairax = gobjects(1,npairs);

        for pr = 1:npairs
            pairax(pr) = nexttile;
            conspair = [grpconds{pr,1} grpconds{pr,2}]; %this panel's conditions

            %fewest trials of any condition in this panel
            nuse = numel(trcond);
            for con = conspair
                if sum(trcond == con) > 0
                    nuse = min(nuse,sum(trcond == con));
                end
            end

            hold on
            for con = conspair
                idxcon = find(trcond == con);
                if isempty(idxcon) == true
                    continue %this animal doesn't have this condition
                end
                pick = idxcon(randperm(numel(idxcon),nuse)); %nuse of this condition's trials
                meancon = reshape(mean(pcscore{m,csdrms_s}(pick,:,:),1),nb,pca_dims); %[plot bins x pca_dims]
                lims = [min(lims(1),min(meancon(:))) max(lims(2),max(meancon(:)))];
                if pca_dims == 3
                    plot3(meancon(:,1),meancon(:,2),meancon(:,3),'Color',condcol(con,:),'LineWidth',1.5,'DisplayName',condis{con})
                    plot3(meancon(onset_bin,1),meancon(onset_bin,2),meancon(onset_bin,3),'o','Color',condcol(con,:),'MarkerFaceColor',condcol(con,:),'HandleVisibility','off')
                else
                    plot(meancon(:,1),meancon(:,2),'Color',condcol(con,:),'LineWidth',1.5,'DisplayName',condis{con})
                    plot(meancon(onset_bin,1),meancon(onset_bin,2),'o','Color',condcol(con,:),'MarkerFaceColor',condcol(con,:),'HandleVisibility','off')
                end
            end
            hold off
            legend('Location','best','Interpreter','none')
            title([char(pairnames(pr)) ', ' num2str(nuse) ' trials each'],'Interpreter','none')
            xlabel(['PC1 (' num2str(pcexpl{m,csdrms_s}(1),'%.1f') '%)'])
            ylabel(['PC2 (' num2str(pcexpl{m,csdrms_s}(2),'%.1f') '%)'])
            if pca_dims == 3
                zlabel(['PC3 (' num2str(pcexpl{m,csdrms_s}(3),'%.1f') '%)'])
                view(3)
                grid on
            end
        end

        %same axes on every panel
        for pr = 1:npairs
            xlim(pairax(pr),lims)
            ylim(pairax(pr),lims)
            if pca_dims == 3
                zlim(pairax(pr),lims)
            end
        end
        sgtitle([char(mice(m)) ' ' char(region) ' CSD RMS: mean of each condition, dot = stim onset'],'Interpreter','none')

        if plotsave == 1
            fname = sprintf([char(mice(m)) '-' 'CSDRMS_byIntensity']);
            savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
        end
    end
end


%% Data Analysis: participation ratio of each condition, split into signal and noise

%PR = (sum eig)^2 / sum(eig^2), about how many dimensions the variance is spread over. This uses all
%the channels, not just the plotted PCs. Three covariances for each condition, each window:
%   Total  - all trials end to end
%   Noise  - each trial minus its condition's mean, so what differs from trial to trial
%   Signal - the part of the mean that repeats: split the trials into two random halves, and take the
%            covariance BETWEEN the two halves' means. Noise is different in each half, so it averages
%            away here. Using the plain trial mean instead would count leftover noise as signal
%SignalFrac = how much of the variance is signal. When it's near 0 (like the baseline should be) there's
%hardly any signal to measure, so ignore PR_signal for that row.
%within an animal every condition uses the same number of trials, because PR depends on trial count
if run_pr == 1
    prmeas = ["Total","Signal","Noise"];
    prmat = nan(nmice,ncond,nsig,nwin,3); %[animals x conditions x signals x windows x Total/Signal/Noise]
    prfrac = nan(nmice,ncond,nsig,nwin); %[animals x conditions x signals x windows] SignalFrac
    prres = table();

    tic
    for m = 1:nmice
        for s = 1:nsig
            %which conditions does this animal have, and what's the fewest trials?
            havecond = [];
            nuse = inf;
            for con = 1:ncond
                if isempty(stimdat{m,con,s}) == false
                    havecond(end+1) = con; %#ok<SAGROW>
                    nuse = min(nuse,size(stimdat{m,con,s},1));
                end
            end
            if isempty(havecond) == true
                continue %no files for this animal
            end
            disp([char(signames(s)) ': PR for ' char(mice(m))])

            %SD of each channel over all this animal's conditions, so a gain change isn't normalised away
            pooled = cat(1, stimdat{m,havecond,s}); %[all trials x nbins x channels]
            chgood = all(isfinite(pooled),[1 2]); %channels with no CSD sit this signal out
            nch = nnz(chgood);
            sdch = ones(1,1,nch);
            if znorm == 1
                sdch = reshape(std(reshape(pooled(:,:,chgood),[],nch),0,1),1,1,nch);
            end
            clear pooled

            for w = 1:nwin
                inwin = time_bins >= win_ms(w,1) & time_bins < win_ms(w,2);
                nb = nnz(inwin);
                for con = havecond
                    X = double(stimdat{m,con,s}(:,inwin,chgood)) ./ sdch; %[trials x window bins x channels]
                    ntr = size(X,1);

                    %add up the covariances over random picks and splits of trials, then take the eigenvalues once
                    Cs = zeros(nch,nch,3); %Total, Signal, Noise
                    for it = 1:pr_nsub
                        pick = randperm(ntr,nuse); %nuse of this condition's trials, in a random order
                        half = floor(nuse/2);
                        meanA = reshape(mean(X(pick(1:half),:,:),1),nb,nch); %[window bins x channels]
                        meanB = reshape(mean(X(pick((half+1):nuse),:,:),1),nb,nch);
                        meanA = meanA - mean(meanA,1); %take out the mean over time
                        meanB = meanB - mean(meanB,1);
                        Cs(:,:,2) = Cs(:,:,2) + (meanA'*meanB + meanB'*meanA)/(2*(nb-1)); %signal: between the two halves

                        Xpick = X(pick,:,:);
                        Cs(:,:,1) = Cs(:,:,1) + cov(reshape(Xpick,[],nch)); %total
                        Cs(:,:,3) = Cs(:,:,3) + cov(reshape(Xpick - mean(Xpick,1),[],nch)); %noise: each trial minus the mean
                    end
                    Cs = Cs/pr_nsub;

                    %PR from each covariance. Signal can have small negative eigenvalues from sampling, count those as 0
                    for mm = 1:3
                        lam = max(eig((Cs(:,:,mm)+Cs(:,:,mm)')/2),0);
                        prmat(m,con,s,w,mm) = sum(lam)^2/sum(lam.^2);
                    end
                    sigfrac = max(trace(Cs(:,:,2)),0)/trace(Cs(:,:,1)); %how much of the variance is signal
                    prfrac(m,con,s,w) = sigfrac;

                    if con <= numel(condis_L) %which group is this condition in?
                        grpname = "L";
                    else
                        grpname = "LW";
                    end
                    prres = [prres; table(mice(m),signames(s),win_names(w),string(condis{con}),grpname,nch,nuse, ...
                        prmat(m,con,s,w,1),prmat(m,con,s,w,2),prmat(m,con,s,w,3),sigfrac, ...
                        'VariableNames',{'Animal_Name','DataType','Window','Condition_Name','Group','nChans','nTrials_used', ...
                        'PR_total','PR_signal','PR_noise','SignalFrac'})]; %#ok<AGROW>
                end
            end
        end
    end
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: participation ratio across light intensity

%one figure per window, rows = Total/Signal/Noise, columns = signals. Thin lines = animals, thick = their mean
if run_pr == 1
    groupconds = {1:numel(condis_L), numel(condis_L) + (1:numel(condis_LW))}; %L conditions, LW conditions
    groupcols = {col_L, col_LW};
    for w = 1:nwin
        figure()
        set(gcf, 'Position', get(0, 'Screensize'));
        tiledlayout(3,nsig)
        sgtitle([char(region) ' ' char(win_names(w)) ', ' num2str(win_ms(w,1)) ' to ' num2str(win_ms(w,2)) ' ms: participation ratio, L = black, LW = red'],'Interpreter','none')
        for mm = 1:3
            for s = 1:nsig
                nexttile
                hold on
                for g = 1:2
                    cons = groupconds{g};
                    prvals = prmat(:,cons,s,w,mm); %[animals x conditions in this group]
                    for m = 1:nmice
                        plot(1:numel(cons),prvals(m,:),'-','Color',[groupcols{g} 0.35])
                    end
                    plot(1:numel(cons),mean(prvals,1,'omitnan'),'-o','Color',groupcols{g},'MarkerFaceColor',groupcols{g},'LineWidth',2)
                end
                hold off
                xlim([0.5 max(numel(condis_L),numel(condis_LW))+0.5])
                xticks(1:numel(condis_L))
                xticklabels(condis_L)
                set(gca,'TickLabelInterpreter','none')
                xlabel('Light intensity (L names shown)')
                ylabel('PR (dims)')
                title([char(signames(s)) ', ' char(prmeas(mm))],'Interpreter','none')
            end
        end

        if plotsave == 1
            fname = sprintf([char(region) '-' char(win_names(w)) '-' 'PR_signal_noise']);
            savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
        end
    end
end


%% Data Analysis: signal PR vs noise PR

%is the repeatable response lower-dimensional than the trial-to-trial variability? Signal and noise come from
%the same trials, so there's nothing to shuffle here. These tests use animals as the n, and only the
%conditions with SignalFrac >= pr_minfrac, because PR_signal means little without signal.
%signal PR is measured from the trial means over time, so it's capped by how many time bins the window has
%(fewer bins = lower cap), and noise PR isn't. So compare windows with care, especially the short P1
if run_pr == 1
    stats_pr_signoise = table(); %one row per signal x window x condition
    stats_pr_lme = table(); %one row per signal x window
    pr_lme_models = cell(nsig,nwin);

    for s = 1:nsig
        for w = 1:nwin
            %one row per animal x condition x Signal/Noise, for the mixed model
            lmeAnimal = strings(0,1);
            lmeMeasure = strings(0,1);
            lmeGroup = strings(0,1);
            lmePR = zeros(0,1);

            for con = 1:ncond
                keepm = find(prfrac(:,con,s,w) >= pr_minfrac); %animals with enough signal in this condition
                if con <= numel(condis_L) %which group is this condition in?
                    grpname = "L";
                else
                    grpname = "LW";
                end
                nhave = sum(isnan(prmat(:,con,s,w,2)) == false);
                if isempty(keepm) == true
                    stats_pr_signoise = [stats_pr_signoise; table(signames(s),win_names(w),string(condis{con}),nhave,0,NaN,NaN,NaN,NaN, ...
                        'VariableNames',{'DataType','Window','Condition_Name','nAnimals','nEnoughSignal','PR_signal_med','PR_noise_med','nSignalLower','signrank_p'})]; %#ok<AGROW>
                    continue
                end
                sigvals = prmat(keepm,con,s,w,2);
                noisevals = prmat(keepm,con,s,w,3);
                stats_pr_signoise = [stats_pr_signoise; table(signames(s),win_names(w),string(condis{con}),nhave,numel(keepm),median(sigvals),median(noisevals),sum(sigvals < noisevals),signrank(sigvals,noisevals), ...
                    'VariableNames',{'DataType','Window','Condition_Name','nAnimals','nEnoughSignal','PR_signal_med','PR_noise_med','nSignalLower','signrank_p'})]; %#ok<AGROW>

                for j = 1:numel(keepm)
                    lmeAnimal = [lmeAnimal; mice(keepm(j)); mice(keepm(j))]; %#ok<AGROW>
                    lmeMeasure = [lmeMeasure; "Signal"; "Noise"]; %#ok<AGROW>
                    lmeGroup = [lmeGroup; grpname; grpname]; %#ok<AGROW>
                    lmePR = [lmePR; log(sigvals(j)); log(noisevals(j))]; %#ok<AGROW>
                end
            end

            %mixed model: is signal PR lower than noise PR overall, with each animal getting its own difference?
            %needs 2+ animals. Log, so the effect is a ratio (signal/noise)
            nanimals_lme = numel(unique(lmeAnimal));
            if nanimals_lme < 2
                continue
            end
            tbl = table(categorical(lmeAnimal),categorical(lmeMeasure),categorical(lmeGroup),lmePR,'VariableNames',{'Animal','Measure','Group','logPR'});
            if numel(unique(lmeGroup)) > 1
                pr_lme_models{s,w} = fitlme(tbl,'logPR ~ Measure + Group + (1 + Measure|Animal)');
            else
                pr_lme_models{s,w} = fitlme(tbl,'logPR ~ Measure + (1 + Measure|Animal)');
            end
            [~,~,fx] = fixedEffects(pr_lme_models{s,w},'DFMethod','satterthwaite');
            rowM = strcmp(fx.Name,'Measure_Signal'); %signal, relative to noise
            stats_pr_lme = [stats_pr_lme; table(signames(s),win_names(w),nanimals_lme,numel(lmePR)/2,exp(fx.Estimate(rowM)),exp(fx.Lower(rowM)),exp(fx.Upper(rowM)),fx.DF(rowM),fx.pValue(rowM), ...
                'VariableNames',{'DataType','Window','nAnimals','nConditionRows','SignalOverNoise','CI_lo','CI_hi','DF','p'})]; %#ok<AGROW>
        end
    end
end


%% Data Analysis: L vs LW shuffle test on signal PR, noise PR, and their ratio

%same idea as the path length test: shuffle a pair's trials between L and LW and see how often the difference is
%as big as the real one. Equal trials in L and LW, using all of the pair's matched trials. Each shuffle uses one
%random half split, and so does each of the pr_nsub real picks, so these numbers are a bit different from prres
%(which averages the covariances first). The real difference is averaged over picks, so the test leans cautious.
%q = Benjamini-Hochberg within a signal x window x measure. "All animals" rows only speak for these animals.
%shuffling trials also counts drift between files as a difference: check the Baseline rows
if run_pr == 1 && pr_test == 1
    prtestnames = ["Signal","Noise","Signal/Noise"];
    stats_pr = table();

    tic
    for s = 1:nsig
        nullpairall = zeros(pr_nperm,3,npairs,nwin); %[shuffles x measures x pairs x windows] LW - L summed over animals
        realpairall = zeros(3,npairs,nwin);
        valpairall = zeros(3,2,npairs,nwin); %[measures x L/LW x pairs x windows] values summed over animals
        npairanimals = zeros(npairs,nwin); %how many animals have each pair

        for m = 1:nmice
            %which conditions does this animal have?
            havecond = [];
            for con = 1:ncond
                if isempty(stimdat{m,con,s}) == false
                    havecond(end+1) = con; %#ok<SAGROW>
                end
            end
            if isempty(havecond) == true
                continue %no files for this animal
            end
            disp([char(signames(s)) ': PR shuffle test for ' char(mice(m))])

            %SD of each channel over all this animal's conditions, same as the PR section
            pooled = cat(1, stimdat{m,havecond,s}); %[all trials x nbins x channels]
            chgood = all(isfinite(pooled),[1 2]); %channels with no CSD sit this signal out
            nch = nnz(chgood);
            sdch = ones(1,1,nch);
            if znorm == 1
                sdch = reshape(std(reshape(pooled(:,:,chgood),[],nch),0,1),1,1,nch);
            end
            clear pooled

            for w = 1:nwin
                inwin = time_bins >= win_ms(w,1) & time_bins < win_ms(w,2);
                nb = nnz(inwin);
                nullanimal = zeros(pr_nperm,3); %LW - L summed over this animal's intensities
                realanimal = zeros(1,3);

                for pr = 1:npairs
                    datL = cat(1, stimdat{m,grpconds{pr,1},s}); %this pair's L trials
                    datLW = cat(1, stimdat{m,grpconds{pr,2},s}); %this pair's LW trials
                    if isempty(datL) == true || isempty(datLW) == true
                        continue %this animal doesn't have both
                    end
                    nL = size(datL,1);
                    nLW = size(datLW,1);
                    nuse = min(nL,nLW); %both get this many trials
                    half = floor(nuse/2);
                    ntot = nL + nLW;
                    X = double(cat(1, datL(:,inwin,chgood), datLW(:,inwin,chgood))) ./ sdch; %[L then LW trials x window bins x channels]
                    clear datL datLW

                    %each trial's x'*x, so the noise covariance of any set of trials is a quick sum
                    trxx = zeros(nch,nch,ntot);
                    for tr = 1:ntot
                        x = reshape(X(tr,:,:),nb,nch); %[window bins x channels]
                        trxx(:,:,tr) = x'*x;
                    end

                    %first the real labels, then the shuffles
                    nits = pr_nsub + pr_nperm;
                    prvals = zeros(nits,3,2); %[goes x Signal/Noise/ratio x L/LW]
                    for it = 1:nits
                        if it <= pr_nsub
                            pickgrp = {randperm(nL,nuse), nL + randperm(nLW,nuse)}; %real labels, equal trials
                        else
                            shuf = randperm(ntot); %shuffled labels
                            pickgrp = {shuf(1:nuse), shuf((nuse+1):(2*nuse))};
                        end
                        for g = 1:2
                            idx = pickgrp{g};
                            meanA = reshape(mean(X(idx(1:half),:,:),1),nb,nch); %[window bins x channels]
                            meanB = reshape(mean(X(idx((half+1):nuse),:,:),1),nb,nch);
                            meanAll = reshape(mean(X(idx,:,:),1),nb,nch);
                            meanA = meanA - mean(meanA,1);
                            meanB = meanB - mean(meanB,1);
                            Csig = (meanA'*meanB + meanB'*meanA)/(2*(nb-1)); %signal: between the two halves
                            Cnoise = (sum(trxx(:,:,idx),3) - nuse*(meanAll'*meanAll))/(nuse*nb-1); %noise: each trial minus the mean
                            lamS = max(eig((Csig+Csig')/2),0);
                            lamN = max(eig((Cnoise+Cnoise')/2),0);
                            prS = sum(lamS)^2/sum(lamS.^2);
                            prN = sum(lamN)^2/sum(lamN.^2);
                            prvals(it,:,g) = [prS prN prS/prN];
                        end
                    end
                    clear X trxx

                    %real difference vs the shuffled ones, one table row per measure
                    realvals = reshape(mean(prvals(1:pr_nsub,:,:),1),3,2); %[Signal/Noise/ratio x L/LW]
                    realdiff = mean(prvals(1:pr_nsub,:,2) - prvals(1:pr_nsub,:,1),1); %1 x measures, LW - L
                    nulldiff = prvals((pr_nsub+1):nits,:,2) - prvals((pr_nsub+1):nits,:,1); %[shuffles x measures]
                    for mm = 1:3
                        stats_pr = [stats_pr; table(mice(m),signames(s),win_names(w),pairnames(pr),prtestnames(mm),nuse,realvals(mm,1),realvals(mm,2),realdiff(mm),std(nulldiff(:,mm)), ...
                            (sum(abs(nulldiff(:,mm)) >= abs(realdiff(mm)))+1)/(pr_nperm+1),NaN, ...
                            'VariableNames',{'Animal_Name','DataType','Window','Pair','Measure','nTrials_used','Value_L','Value_LW','Diff','Null_sd','p','q'})]; %#ok<AGROW>
                    end

                    %keep adding up for the "All" rows
                    nullanimal = nullanimal + nulldiff;
                    realanimal = realanimal + realdiff;
                    nullpairall(:,:,pr,w) = nullpairall(:,:,pr,w) + nulldiff;
                    realpairall(:,pr,w) = realpairall(:,pr,w) + realdiff';
                    valpairall(:,:,pr,w) = valpairall(:,:,pr,w) + realvals;
                    npairanimals(pr,w) = npairanimals(pr,w) + 1;
                end

                %this animal, all intensities together
                for mm = 1:3
                    stats_pr = [stats_pr; table(mice(m),signames(s),win_names(w),"All intensities",prtestnames(mm),NaN,NaN,NaN,realanimal(mm),std(nullanimal(:,mm)), ...
                        (sum(abs(nullanimal(:,mm)) >= abs(realanimal(mm)))+1)/(pr_nperm+1),NaN, ...
                        'VariableNames',{'Animal_Name','DataType','Window','Pair','Measure','nTrials_used','Value_L','Value_LW','Diff','Null_sd','p','q'})]; %#ok<AGROW>
                end
            end
        end

        %every animal together: each intensity, then all intensities
        for w = 1:nwin
            for mm = 1:3
                nullwin = zeros(pr_nperm,1);
                realwin = 0;
                for pr = 1:npairs
                    if npairanimals(pr,w) == 0
                        continue %nobody has this pair
                    end
                    nullthis = nullpairall(:,mm,pr,w);
                    realthis = realpairall(mm,pr,w);
                    stats_pr = [stats_pr; table("All animals",signames(s),win_names(w),pairnames(pr),prtestnames(mm),NaN,valpairall(mm,1,pr,w)/npairanimals(pr,w),valpairall(mm,2,pr,w)/npairanimals(pr,w),realthis,std(nullthis), ...
                        (sum(abs(nullthis) >= abs(realthis))+1)/(pr_nperm+1),NaN, ...
                        'VariableNames',{'Animal_Name','DataType','Window','Pair','Measure','nTrials_used','Value_L','Value_LW','Diff','Null_sd','p','q'})]; %#ok<AGROW>
                    nullwin = nullwin + nullthis;
                    realwin = realwin + realthis;
                end
                stats_pr = [stats_pr; table("All animals",signames(s),win_names(w),"All intensities",prtestnames(mm),NaN,NaN,NaN,realwin,std(nullwin), ...
                    (sum(abs(nullwin) >= abs(realwin))+1)/(pr_nperm+1),NaN, ...
                    'VariableNames',{'Animal_Name','DataType','Window','Pair','Measure','nTrials_used','Value_L','Value_LW','Diff','Null_sd','p','q'})]; %#ok<AGROW>
            end
        end
    end

    %correct for doing lots of tests: Benjamini-Hochberg q, separately for single-animal and all-animal rows
    for s = 1:nsig
        for w = 1:nwin
            for mm = 1:3
                for scope = 1:2
                    rowsq = stats_pr.DataType == signames(s) & stats_pr.Window == win_names(w) & stats_pr.Measure == prtestnames(mm) & stats_pr.Pair ~= "All intensities";
                    if scope == 1
                        rowsq = find(rowsq & stats_pr.Animal_Name ~= "All animals");
                    else
                        rowsq = find(rowsq & stats_pr.Animal_Name == "All animals");
                    end
                    if isempty(rowsq) == true
                        continue
                    end
                    [psort,ord] = sort(stats_pr.p(rowsq));
                    nq = numel(psort);
                    qsort = psort.*nq./(1:nq)';
                    for i = (nq-1):-1:1
                        qsort(i) = min(qsort(i),qsort(i+1)); %q can't go down as p goes up
                    end
                    stats_pr.q(rowsq(ord)) = min(qsort,1);
                end
            end
        end
    end
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Visualization: signal PR and noise PR on the same axes

%rows = windows, columns = signals, mean +- SEM over animals. Solid = signal, dashed = noise, black = L, red = LW.
%hollow markers = on average too little signal (SignalFrac < pr_minfrac) to trust that signal PR
if run_pr == 1
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tiledlayout(nwin,nsig)
    sgtitle([char(region) ': signal PR (solid) vs noise PR (dashed), L = black, LW = red, mean +- SEM over animals, hollow = SignalFrac < ' num2str(pr_minfrac)],'Interpreter','none')
    measstyles = {'-','--'}; %signal, noise
    measnames = ["signal","noise"];
    for w = 1:nwin
        for s = 1:nsig
            nexttile
            hold on
            hleg = gobjects(1,4);
            legtxt = strings(1,4);
            nleg = 0;
            for g = 1:2
                cons = groupconds{g};
                xvals = 1:numel(cons);
                fracmean = mean(prfrac(:,cons,s,w),1,'omitnan');
                for j = 1:2
                    prvals = prmat(:,cons,s,w,j+1); %[animals x conditions] Signal, then Noise
                    prmean = mean(prvals,1,'omitnan');
                    prsem = std(prvals,0,1,'omitnan')./sqrt(sum(isnan(prvals) == false,1));
                    nleg = nleg + 1;
                    hleg(nleg) = errorbar(xvals,prmean,prsem,'Color',groupcols{g},'LineStyle',measstyles{j},'Marker','o','MarkerFaceColor','w','LineWidth',1.2);
                    if g == 1
                        legtxt(nleg) = "L " + measnames(j);
                    else
                        legtxt(nleg) = "LW " + measnames(j);
                    end

                    %fill in the markers we trust. Noise PR doesn't need signal, so those are always filled
                    if j == 2
                        trustpts = true(size(xvals));
                    else
                        trustpts = fracmean >= pr_minfrac;
                    end
                    plot(xvals(trustpts),prmean(trustpts),'o','Color',groupcols{g},'MarkerFaceColor',groupcols{g})
                end
            end
            hold off
            xlim([0.5 max(numel(condis_L),numel(condis_LW))+0.5])
            xticks(1:numel(condis_L))
            xticklabels(condis_L)
            set(gca,'TickLabelInterpreter','none')
            ylabel('PR (dims)')
            ylim([0 inf])

            %the mixed model result, if there is one
            titletxt = {[char(signames(s)) ', ' char(win_names(w)) ' (' num2str(win_ms(w,1)) ' to ' num2str(win_ms(w,2)) ' ms)']};
            if height(stats_pr_lme) > 0
                rl = stats_pr_lme(stats_pr_lme.DataType == signames(s) & stats_pr_lme.Window == win_names(w),:);
                if height(rl) == 1
                    if rl.p < 0.001
                        ptxt = 'p < 0.001';
                    else
                        ptxt = ['p = ' num2str(rl.p,'%.3f')];
                    end
                    titletxt{2} = ['signal/noise = ' num2str(rl.SignalOverNoise,'%.2f') ' [' num2str(rl.CI_lo,'%.2f') ' ' num2str(rl.CI_hi,'%.2f') '], ' ptxt];
                end
            end
            title(titletxt,'Interpreter','none','FontSize',9)
            if w == 1 && s == 1
                legend(hleg,cellstr(legtxt),'Location','best')
            end
        end
    end

    if plotsave == 1
        fname = sprintf([char(region) '-' 'PR_signal_vs_noise']);
        savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
    end
end


%% Data Visualization: signal PR / noise PR across light intensity

%under 1 = the repeatable response has fewer dimensions than the trial-to-trial variability.
%thin lines = animals, thick = their mean. * above an intensity = L and LW differ in the all-animal shuffle test (q < 0.05)
if run_pr == 1
    figure()
    set(gcf, 'Position', get(0, 'Screensize'));
    tiledlayout(nwin,nsig)
    sgtitle([char(region) ': signal PR / noise PR, L = black, LW = red, thin = animals, hollow = SignalFrac < ' num2str(pr_minfrac) ', * = L vs LW q < 0.05'],'Interpreter','none')
    for w = 1:nwin
        for s = 1:nsig
            nexttile
            hold on
            allvals = [];
            for g = 1:2
                cons = groupconds{g};
                xvals = 1:numel(cons);
                ratiovals = prmat(:,cons,s,w,2)./prmat(:,cons,s,w,3); %[animals x conditions in this group]
                allvals = [allvals; ratiovals(:)]; %#ok<AGROW>
                for m = 1:nmice
                    plot(xvals,ratiovals(m,:),'-','Color',[groupcols{g} 0.35])
                end
                ratiomean = mean(ratiovals,1,'omitnan');
                plot(xvals,ratiomean,'-o','Color',groupcols{g},'MarkerFaceColor','w','LineWidth',2)
                trustpts = mean(prfrac(:,cons,s,w),1,'omitnan') >= pr_minfrac;
                plot(xvals(trustpts),ratiomean(trustpts),'o','Color',groupcols{g},'MarkerFaceColor',groupcols{g})
            end
            topval = max(allvals);

            %stars from the shuffle test, one per pair
            allintp = NaN;
            if pr_test == 1
                for pr = 1:npairs
                    rq = stats_pr(stats_pr.Animal_Name == "All animals" & stats_pr.DataType == signames(s) & stats_pr.Window == win_names(w) & stats_pr.Measure == "Signal/Noise" & stats_pr.Pair == pairnames(pr),:);
                    if height(rq) == 1 && rq.q < 0.05
                        text(pr,1.05*topval,'*','HorizontalAlignment','center','FontSize',14)
                    end
                end
                rall = stats_pr(stats_pr.Animal_Name == "All animals" & stats_pr.DataType == signames(s) & stats_pr.Window == win_names(w) & stats_pr.Measure == "Signal/Noise" & stats_pr.Pair == "All intensities",:);
                if height(rall) == 1
                    allintp = rall.p;
                end
            end
            hold off
            xlim([0.5 max(numel(condis_L),numel(condis_LW))+0.5])
            xticks(1:numel(condis_L))
            xticklabels(condis_L)
            set(gca,'TickLabelInterpreter','none')
            ylabel('Signal PR / noise PR')
            ylim([0 1.2*topval])
            title({[char(signames(s)) ', ' char(win_names(w))],['L vs LW, all intensities: p = ' num2str(allintp,'%.3f')]},'Interpreter','none','FontSize',9)
        end
    end

    if plotsave == 1
        fname = sprintf([char(region) '-' 'PR_signal_noise_ratio']);
        savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
    end
end



%Now beginning the optional stats section, everything from here to the report needs run_stats = 1

%% Data Analysis (optional): PR, alignment, and angles in L vs LW, each against shuffled labels

if run_stats == 1
    stats_win = table();
    trajdat = cell(nmice,npairs,nsig); %{animal x pair x signal} z-scored data and trial labels, kept for the plots below

    tic
    for s = 1:nsig
        for m = 1:nmice
            for pr = 1:npairs
                datL = cat(1, stimdat{m,grpconds{pr,1},s}); %all L trials of this animal, [trials x nbins x channels]
                datLW = cat(1, stimdat{m,grpconds{pr,2},s}); %all LW trials
                nL = size(datL,1);
                nLW = size(datLW,1);
                nuse = min(nL,nLW); %both conditions get this many trials
                if nuse < min_trials
                    run_notes(end+1) = mice(m) + ", " + pairnames(pr) + ", " + signames(s) + ": no stats, only " + nL + " L and " + nLW + " LW trials"; %#ok<SAGROW>
                    continue
                end
                disp([char(signames(s)) ': Manifold stats for ' char(mice(m)) ', ' char(pairnames(pr))])

                %stick L and LW together
                pooled = cat(1, datL, datLW); %[L then LW trials x nbins x channels]
                clear datL datLW

                %channels with no CSD sit this signal out (noted in the PCA section above)
                pooled = pooled(:,:,all(isfinite(pooled),[1 2]));
                nch = size(pooled,3);

                %z-score each channel by its SD over both conditions
                if znorm == 1
                    sdch = std(reshape(pooled,[],nch),0,1); %1 x channels
                    pooled = pooled ./ reshape(sdch,1,1,nch);
                end
                ntot = nL + nLW;
                trajdat{m,pr,s} = struct('pooled',pooled,'nL',nL,'nuse',nuse);

                %real L vs LW only needs one go if the trial counts already match
                nsub_use = n_sub;
                if nL == nLW
                    nsub_use = 1;
                end
                nits = nsub_use + n_perm; %first the real labels, then the shuffles

                for w = 1:nwin
                    inwin = time_bins >= win_ms(w,1) & time_bins < win_ms(w,2);
                    X = double(pooled(:,inwin,:)); %[trials x window bins x channels]
                    nb = size(X,2);
                    Nsamp = nuse*nb; %samples behind each single-trial covariance

                    %per-trial building blocks, so the single-trial covariance of any set of trials is a quick sum
                    trxx = zeros(nch,nch,ntot); %[channels x channels x trials] x'*x
                    trsum = zeros(ntot,nch); %[trials x channels] sum over time
                    for tr = 1:ntot
                        x = reshape(X(tr,:,:),nb,nch); %[window bins x channels]
                        trxx(:,:,tr) = x'*x;
                        trsum(tr,:) = sum(x,1);
                    end

                    %results of every go, [gos x spaces (x L/LW)]
                    prs = zeros(nits,2,2);
                    totvar = zeros(nits,2,2);
                    aligns = zeros(nits,2);
                    angles = zeros(nits,2);

                    for it = 1:nits
                        if it <= nsub_use
                            pickgrp = {randperm(nL,nuse), nL + randperm(nLW,nuse)}; %real labels, equal trials
                        else
                            shuf = randperm(ntot); %shuffled labels
                            pickgrp = {shuf(1:nuse), shuf((nuse+1):(2*nuse))};
                        end

                        Cs = cell(2,2); %{space x group} channel covariance
                        for g = 1:2
                            idx = pickgrp{g};
                            Cs{1,g} = cov(reshape(mean(X(idx,:,:),1),nb,nch)); %TrialAvg: covariance of the trial mean over time
                            mu = sum(trsum(idx,:),1)/Nsamp; %SingleTrial: covariance over all trials end to end
                            Cs{2,g} = (sum(trxx(:,:,idx),3) - Nsamp*(mu'*mu))/(Nsamp-1);
                        end

                        for sp = 1:2
                            Vk = cell(1,2); %top pc_k PCs of each group
                            lams = cell(1,2); %eigenvalues of each group, biggest first
                            for g = 1:2
                                [V,D] = eig((Cs{sp,g}+Cs{sp,g}')/2);
                                [lam,ord] = sort(diag(D),'descend');
                                lam = max(lam,0); %tiny negatives are rounding error
                                Vk{g} = V(:,ord(1:pc_k));
                                lams{g} = lam;
                                prs(it,sp,g) = sum(lam)^2/sum(lam.^2); %participation ratio
                                totvar(it,sp,g) = sum(lam);
                            end
                            %how much of LW fits in L's top PCs, relative to the best it could do, and vice versa
                            a_LinLW = trace(Vk{2}'*Cs{sp,1}*Vk{2})/sum(lams{1}(1:pc_k));
                            a_LWinL = trace(Vk{1}'*Cs{sp,2}*Vk{1})/sum(lams{2}(1:pc_k));
                            aligns(it,sp) = (a_LinLW + a_LWinL)/2;
                            angles(it,sp) = mean(acosd(min(svd(Vk{1}'*Vk{2}),1))); %mean principal angle, deg
                        end
                    end

                    %real numbers vs the shuffles, one table row per space
                    for sp = 1:2
                        pr_diff_real = mean(prs(1:nsub_use,sp,2) - prs(1:nsub_use,sp,1));
                        pr_diff_null = prs((nsub_use+1):nits,sp,2) - prs((nsub_use+1):nits,sp,1);
                        align_real = mean(aligns(1:nsub_use,sp));
                        align_null = aligns((nsub_use+1):nits,sp);
                        angle_real = mean(angles(1:nsub_use,sp));
                        angle_null = angles((nsub_use+1):nits,sp);

                        row = table(mice(m),signames(s),region,pairnames(pr),win_names(w),spacenames(sp), ...
                            nch,nL,nLW,nuse, ...
                            mean(prs(1:nsub_use,sp,1)),mean(prs(1:nsub_use,sp,2)),pr_diff_real,std(pr_diff_null), ...
                            (sum(abs(pr_diff_null) >= abs(pr_diff_real))+1)/(n_perm+1), ...
                            mean(totvar(1:nsub_use,sp,1)),mean(totvar(1:nsub_use,sp,2)), ...
                            align_real,median(align_null),(sum(align_null <= align_real)+1)/(n_perm+1), ...
                            angle_real,median(angle_null),(sum(angle_null >= angle_real)+1)/(n_perm+1), ...
                            'VariableNames',{'Animal_Name','DataType','Region','Pair','Window','Space', ...
                            'nChans','nTrials_L','nTrials_LW','nTrials_used', ...
                            'PR_L','PR_LW','PR_diff','PR_null_sd','PR_p', ...
                            'TotVar_L','TotVar_LW', ...
                            'Align','Align_null','Align_p', ...
                            'Angle','Angle_null','Angle_p'});
                        stats_win = [stats_win; row]; %#ok<AGROW>
                    end
                end
                clear pooled X trxx trsum
            end
        end
    end
    elapsed_time = toc;
    fprintf('execution time: %.2f seconds\n', elapsed_time);
end


%% Data Analysis (optional): across animals, paired L vs LW

if run_stats == 1
    stats_paired = table();
    for s = 1:nsig
        for pr = 1:npairs
            for w = 1:nwin
                for sp = 1:2
                    r = stats_win(stats_win.DataType == signames(s) & stats_win.Pair == pairnames(pr) & stats_win.Window == win_names(w) & stats_win.Space == spacenames(sp),:);
                    if height(r) == 0
                        continue
                    end
                    row = table(signames(s),pairnames(pr),win_names(w),spacenames(sp),height(r), ...
                        median(r.PR_L),median(r.PR_LW),sum(r.PR_LW < r.PR_L),signrank(r.PR_L,r.PR_LW), ...
                        median(r.Align),median(r.Align_null), ...
                        'VariableNames',{'DataType','Pair','Window','Space','nAnimals','PR_L_med','PR_LW_med','nAnimals_LWlower','PR_signrank_p','Align_med','Align_null_med'});
                    stats_paired = [stats_paired; row]; %#ok<AGROW>
                end
            end
        end
    end
end


%% Data Visualization (optional): L and LW trial means in shared PC space

if run_stats == 1 && plot_traj == 1
    trajL = cell(nmice,npairs,nsig); %{animal x pair x signal} [plot bins x pc_k] shared-PC scores
    trajLW = cell(nmice,npairs,nsig);
    jointexpl = cell(nmice,npairs,nsig);

    disp('Manifold: Making some more pictures for you')
    for s = 1:nsig
        for m = 1:nmice
            for pr = 1:npairs
                if isempty(trajdat{m,pr,s}) == true
                    continue %skipped above, too few trials
                end
                pooled = trajdat{m,pr,s}.pooled(:,inplot,:);
                nL = trajdat{m,pr,s}.nL;
                nuse = trajdat{m,pr,s}.nuse;
                nb = size(pooled,2);
                nch = size(pooled,3);

                %same number of trials in both trial means, so one isn't noisier than the other
                avgL = reshape(mean(pooled(randperm(nL,nuse),:,:),1),nb,nch); %[plot bins x channels]
                avgLW = reshape(mean(pooled(nL + randperm(size(pooled,1)-nL,nuse),:,:),1),nb,nch);

                %fit the PCs to both conditions together, so both live in the same axes
                [~,score,latent] = pca([avgL; avgLW]); %score is [2*bins x PCs]
                jointexpl{m,pr,s} = 100*latent/sum(latent); %percent variance per PC
                trajL{m,pr,s} = score(1:nb,1:pc_k);
                trajLW{m,pr,s} = score((nb+1):(2*nb),1:pc_k);

                figure()
                set(gcf, 'Position', get(0, 'Screensize'));
                tiledlayout(2,3)
                sgtitle([char(mice(m)) ' ' char(region) ' ' char(signames(s)) ': ' char(pairnames(pr)) ', trial means, ' num2str(nuse) ' trials each'],'Interpreter','none')

                %the manifold itself, first 3 shared PCs
                nexttile
                plot3(trajL{m,pr,s}(:,1),trajL{m,pr,s}(:,2),trajL{m,pr,s}(:,3),'Color',col_L)
                hold on
                plot3(trajLW{m,pr,s}(:,1),trajLW{m,pr,s}(:,2),trajLW{m,pr,s}(:,3),'Color',col_LW)
                plot3(trajL{m,pr,s}(onset_bin,1),trajL{m,pr,s}(onset_bin,2),trajL{m,pr,s}(onset_bin,3),'o','Color',col_L,'MarkerFaceColor',col_L)
                plot3(trajLW{m,pr,s}(onset_bin,1),trajLW{m,pr,s}(onset_bin,2),trajLW{m,pr,s}(onset_bin,3),'o','Color',col_LW,'MarkerFaceColor',col_LW)
                hold off
                grid on
                view(3)
                xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
                legend({'L','LW','stim onset (L)','stim onset (LW)'},'Location','best')
                title('Shared PC space')

                %first 3 shared PCs over time
                for pc = 1:3
                    nexttile
                    plot(tplot,trajL{m,pr,s}(:,pc),'Color',col_L)
                    hold on
                    plot(tplot,trajLW{m,pr,s}(:,pc),'Color',col_LW)
                    hold off
                    xline(0)
                    xlim(plot_ms)
                    xlabel('ms re: stim')
                    title(['PC' num2str(pc) ' (' num2str(jointexpl{m,pr,s}(pc),'%.1f') '% var)'])
                end

                %eigenvalue spectra of each condition on its own, plot stretch only
                lamL = sort(eig(cov(avgL)),'descend');
                lamLW = sort(eig(cov(avgLW)),'descend');
                lamL = max(lamL,0);
                lamLW = max(lamLW,0);
                prL_plot = sum(lamL)^2/sum(lamL.^2); %PR of what's in this plot
                prLW_plot = sum(lamLW)^2/sum(lamLW.^2);
                npcs_show = min(15,nch);
                nexttile
                plot(1:npcs_show,100*cumsum(lamL(1:npcs_show))/sum(lamL),'-o','Color',col_L)
                hold on
                plot(1:npcs_show,100*cumsum(lamLW(1:npcs_show))/sum(lamLW),'-o','Color',col_LW)
                hold off
                ylim([0 100])
                xlabel('PCs (each condition on its own)')
                ylabel('Cumulative % variance')
                legend({['L, PR = ' num2str(prL_plot,'%.1f')],['LW, PR = ' num2str(prLW_plot,'%.1f')]},'Location','southeast')
                title('Dimensionality, trial mean')

                %how far apart are the two trajectories at each time point?
                nexttile
                plot(tplot,sqrt(sum((trajL{m,pr,s} - trajLW{m,pr,s}).^2,2)),'k')
                xline(0)
                xlim(plot_ms)
                xlabel('ms re: stim')
                ylabel(['L-LW distance, top ' num2str(pc_k) ' PCs'])
                title('Trajectory separation')

                if plotsave == 1
                    fname = sprintf([char(mice(m)) '-' char(signames(s)) '-' char(pairnames(pr)) '-' 'manifold_trajectories']);
                    savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
                end
            end
        end
    end
end


%% Data Visualization (optional): paired L vs LW across animals

if run_stats == 1 && plot_summary == 1 && height(stats_win) > 0
    %what each figure shows: {column A, column B, x labels, y label, file name part}
    sumfigs = {'PR_L','PR_LW',{'L','LW'},'Participation ratio','PR'; ...
               'Align_null','Align',{'Shuffled','Real'},'Alignment','Alignment'};

    for s = 1:nsig
        for pr = 1:npairs
            for f = 1:size(sumfigs,1)
                figure()
                set(gcf, 'Position', get(0, 'Screensize'));
                tiledlayout(2,nwin)
                sgtitle([char(region) ' ' char(signames(s)) ': ' char(pairnames(pr)) ', ' sumfigs{f,4}],'Interpreter','none')
                for sp = 1:2
                    for w = 1:nwin
                        nexttile
                        r = stats_win(stats_win.DataType == signames(s) & stats_win.Pair == pairnames(pr) & stats_win.Window == win_names(w) & stats_win.Space == spacenames(sp),:);
                        hold on
                        for i = 1:height(r) %one line per animal
                            plot([1 2],[r.(sumfigs{f,1})(i) r.(sumfigs{f,2})(i)],'-o','Color',[0.5 0.5 0.5])
                        end
                        plot([1 2],[median(r.(sumfigs{f,1})) median(r.(sumfigs{f,2}))],'-','Color',col_LW,'LineWidth',2) %median across animals
                        hold off
                        xlim([0.5 2.5])
                        xticks([1 2])
                        xticklabels(sumfigs{f,3})
                        ylabel(sumfigs{f,4})
                        title([char(win_names(w)) ', ' char(spacenames(sp))])
                    end
                end

                if plotsave == 1
                    fname = sprintf([char(signames(s)) '-' char(pairnames(pr)) '-' 'manifold_' sumfigs{f,5}]);
                    savefig(gcf,[save_directory '\Figures_Manifold\' fname],'compact')
                end
            end
        end
    end
end


%% Save the stats as a csv for R (optional)

if run_stats == 1 && save_csv == 1
    [csvname,csvpath] = uiputfile(fullfile(csv_directory,'ManifoldPCA_LvsLW.csv'),'Name the csv');
    if ischar(csvname)
        writetable(stats_win,fullfile(csvpath,csvname));
        disp('saved')
    else
        run_notes(end+1) = "csv not saved, the name box was cancelled";
    end
end


%% Post-run report

disp('---------- Manifold report ----------')
disp(filesused)
disp(pathres)
disp(stats_path)
disp(stats_path_paired)
disp(stats_path_lme)
disp('Path test: shuffles trials, so drift between files also counts. Check it with path_ms before the stim. p is uncorrected, q is corrected within a signal.')
if run_pr == 1
    disp(prres)
    disp(stats_pr_signoise)
    disp(stats_pr_lme)
    if pr_test == 1
        disp(stats_pr)
    end
    disp('PR: only trust PR_signal where SignalFrac is well above 0. PR depends on channel count, so compare within an animal.')
end
if run_stats == 1
    disp(stats_win(:,{'Animal_Name','DataType','Pair','Window','Space','nChans','nTrials_used','PR_L','PR_LW','PR_diff','PR_p','Align','Align_null','Align_p','Angle','Angle_null','Angle_p'}))
    disp(stats_paired)
    disp('P values are uncorrected across signals, windows and measures.')
end
for i = 1:numel(run_notes)
    disp(run_notes(i))
end
