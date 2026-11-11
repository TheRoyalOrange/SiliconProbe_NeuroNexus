%% OscillationEvents_StimPrePost.m
%
% Description: Runs oshkosh_detect_funct on every recording you pick for one
%   animal, sorts the probe events (prbEvent) by whether a stimulus was
%   present, and measures the participation ratio (PR) across each probe's
%   channels, the same way NNx_ManifoldPCA_LvsLW.m does: on the Hilbert
%   amplitude envelope of the LFP in each band of hbert_bands (default 4-12,
%   12-30 and 30-80 Hz), binned at bin_ms and z-scored per channel.
%   PR = (sum of eigenvalues)^2 / (sum of squared eigenvalues) of the channel
%   covariance matrix. It's about how many independent dimensions the activity
%   uses: 1 = every channel doing the same thing, number of channels = every
%   channel doing its own thing.
%   PR is measured THREE WAYS so they can be compared:
%     per-window - one PR per event / stim window, then averaged. A short window
%       of a narrow-band envelope only has a few independent time points, and PR
%       can't go above that, so this is capped in the slow bands
%     pooled - every window of a group (mean-centred first) stacked into ONE
%       covariance, then one PR. Many more time points, so no cap, but no
%       per-window values, and loud windows count for more (their covariance is bigger)
%     equal-weight pooled - same as pooled, but each window's covariance is scaled
%       to the same total variance (trace = 1) first, so every window counts the
%       same, like per-window does, while keeping pooling's many time points
%   for: each event, pre vs post stim windows (pooled gets a swap-label
%   permutation test, per-window gets signrank), and a sliding window over time
%   (aligned to stim onset, and to event start). Also checks how many events go
%   with each stim, draws the events as lines through PC space (faded = early,
%   full colour = late), draws each recording's stims as lines through PC space
%   lined up on stim onset (grey = before, red = after), and opens the
%   scrollable event viewer
%   (oscillation_plotEvents). Nothing is saved to disk, everything stays in the
%   workspace.
%
% Inputs:
%   animal (string) - animal ID, e.g. '20260821-p10'. Recording folder names
%     must contain <animal>-<condition> (oshkosh_detect_funct errors otherwise)
%   save_directory (string) - folder holding <animal>-ProbeInfo.mat
%   fq_range (double, 1x2, Hz) - bandpass used for oscillation detection. The
%     snippets PR is measured from are filtered in this band, so hbert_bands
%     have to sit inside it
%   probes (int vector) - indices into ProbeInfo.ProbeIds/ProbeMaps to analyze
%   reorder (0/1) - pre-06/2024 channel remap flag for oshkosh_detect_funct
%   nostim_conds (cell of strings) - recordings whose condition name contains
%     one of these have their stim_times and during_stim flags wiped
%   stim_reach (double, ms) - an event starting this long after a stim still
%     counts as "with stim". 3000 matches the viewer's grey shading
%   train_gap (double, ms) - TTL-high samples closer than this are one stim
%   pre_win, post_win (double, 1x2, ms) - windows relative to stim onset
%   hbert_bands (double, bands x 2, Hz) - one Hilbert envelope per row
%   bin_ms (double, ms) - bin width of the envelopes
%   znorm (0/1) - z-score each channel by its SD over all event bins
%   edge_ms (double, ms) - a window is only used if its snippet covers it with
%     this much to spare on each side (filter and Hilbert edge wobble)
%   n_perm (double) - label swaps for the pooled pre vs post test
%   rng_seed (double) - same swaps every run
%   slide_ms, slide_step_ms (double, ms) - sliding window length and step
%   slide_stim_ms, slide_event_ms (double, 1x2, ms) - stretch the sliding
%     window covers, re: stim onset and re: event start
%   plot_pca, pca_dims, ev_plot_ms, fade_start, trial_alpha - PCA plot settings
%   plot_stimalign (0/1) - draw the stim-aligned PC plot, one figure per recording
%   align_pcs (int vector, 1x2 or 1x3) - which PCs the stim-aligned plot shows
%     (e.g. [2 3 4]). 2 = flat plot, 3 = 3D. Doesn't use pca_dims
%   align_win (double, 1x2, ms) - stretch around stim onset each stim's line covers
%   align_alpha, align_width (double) - see-through-ness and thickness of each stim's line
%   align_shift (0/1) - 1 = move each line so its onset point sits at 0
%   align_pcfit (0/1) - 0 = one PCA per probe over all recordings, 1 = one per recording
%   view_which (0/1/2) - which events the scrollable viewers show
%   view_pr (0/1) - open the bout viewer with PR over time and a 3-PC manifold
%   ProbeInfo (struct) - loaded from <save_directory>\<animal>-ProbeInfo.mat
%     (see oshkosh_detect_funct.m and oscillation_plotEvents.m for fields used)
%   data_directory (cell of strings, picked via uipickfiles, made Nfiles x 1) -
%     raw OpenEphys recording folders
%
% Outputs (workspace only):
%   oscillations, snippets, prbsnippets (structs) - from oshkosh_detect_funct.
%     oscillations has its stim_times and during_stim fields wiped for
%     nostim_conds recordings
%   fileReport (table, Nfiles rows) - condition, TTL-high samples oshkosh found
%     (before any wiping) and stim onsets kept, per recording. Check this after
%     the run: n_ttl_samples above 0 on a spon recording means TTL noise was
%     read as stim. n_onsets should match the stims you actually gave
%   stim_file, stim_onset, stim_condition (Nstim x 1) - recording index, onset
%     (1kHz samples, local to that recording) and condition of every stim
%   signames (string, 1 x bands) - "LFP_hbert_<lo>to<hi>", one per hbert_bands row
%   bandtags (string, 1 x bands) - "<lo>to<hi>", the band part of table column names
%   envz (cell, 1 x probes) - {probe}{band}{event} (single, bins x channels)
%     Hilbert envelope of the whole snippet (event +/- 1 s), binned, z-scored
%     (z-units if znorm = 1, uV otherwise)
%   binstart (cell, 1 x probes) - {probe}{event} (double, bins x 1) recording
%     sample (1 kHz, local to that recording) where each bin of envz starts
%   evmask (cell, 1 x probes) - {probe}{event} (logical, bins x 1) the bins
%     inside the event itself
%   sdch (cell, 1 x probes) - (double, bands x channels) SD used for the z-score
%   eventTable (cell, 1 x length(probes)) - each a table, one row per prbEvent:
%     file, condition, start, stop (samples), dur_ms, peak_rms (uV),
%     during_stim (TTL overlap only), with_stim (TTL overlap or within
%     stim_reach), n_bins (bins inside the event), pr_<band> (per-window PR
%     over the event, 1 to Nch, one column per band)
%   stimSortSummary (cell, 1 x length(probes)) - table, one row per stim group
%     (no stim / with stim): n_events, medians of dur_ms and peak_rms, and per
%     band med_pr_<band> (median per-event PR), poolpr_<band> (pooled PR) and
%     pooleqpr_<band> (equal-weight pooled PR)
%   stim_ev (cell, 1 x probes) - (double, Nstim x 1) the event that stim's
%     onset lands inside, if its snippet covers the pre AND post window with
%     edge_ms to spare, NaN if none
%   pr_pre, pr_post (cell, 1 x length(probes)) - each (double, Nstim x bands):
%     per-window PR in that window for each stim, NaN if not covered
%   prepostSummary (cell, 1 x length(probes)) - table, one row per stim
%     condition: n_stim, n_used (stims with both windows covered), and per band:
%     prwin_pre_<band>, prwin_post_<band> (mean per-window PR), signrank_p_<band>
%     (post vs pre, per-window), prpool_pre_<band>, prpool_post_<band> (pooled
%     PR), perm_p_<band> (swap-label test on the pooled post - pre, two-sided),
%     prpooleq_pre_<band>, prpooleq_post_<band>, permeq_p_<band> (the same for
%     equal-weight pooled),
%     n_indep_<band> (rough number of independent time points behind each pooled
%     PR: n_used x window length x band width. Should be well above nChans)
%   slide_t_stim, slide_t_ev (double, 1 x steps) - sliding window centres, ms re:
%     stim onset / re: event start
%   slidePR_stim_each (cell, 1 x probes) - (double, Nstim x steps x bands)
%     per-window PR of each stim's sliding windows, NaN if not covered
%   slidePR_stim_pool (cell, 1 x probes) - (double, stim conditions x steps x
%     bands) pooled PR across that condition's stims at each step
%   slidePR_stim_pooleq (cell, 1 x probes) - same, equal-weight pooled
%   slide_n_stim (cell, 1 x probes) - (double, stim conditions x steps) stims behind each step
%   slidePR_ev_each (cell, 1 x probes) - (double, events x steps x bands)
%   slidePR_ev_pool (cell, 1 x probes) - (double, 2 groups x steps x bands),
%     group 1 = no stim, 2 = with stim
%   slidePR_ev_pooleq (cell, 1 x probes) - same, equal-weight pooled
%   slide_n_ev (cell, 1 x probes) - (double, 2 x steps) events behind each step
%   stimEventReport (cell, 1 x probes) - table, one row per recording with
%     stims: n_stims, stims with 0 / 1 / 2+ linked events (event starts within
%     stim_reach after the onset, or spans it), events linked, and events
%     linked to more than one stim
%   alignTraj (cell, probes x recordings) - {band} (double, stims x bins x
%     numel(align_pcs)) each stim's line through the align_pcs PCs over align_win, shifted to 0
%     at onset if align_shift = 1. Empty for recordings with no stim inside an event
%   alignStim (cell, probes x recordings) - (double, stims x 1) which stims
%     (index into stim_onset) the rows of alignTraj are
%   run_notes (string) - problems hit along the way, shown in the report at the end
%   Figures - one "all events" figure per recording (from oshkosh), and per
%     probe: events, pre vs post, sliding PR and PCA trajectory figures, one
%     stim-aligned PCA figure per recording with stims inside events, plus
%     two scrollable viewers: oscillation_plotEvents, and (view_pr = 1) the bout
%     viewer: every probe's bouts in one time-ordered list. Per bout, over
%     exactly that bout's snippet: sliding per-window PR in every band, each
%     probe on its own channels (top, bout's probe solid, other probes dashed),
%     and per shown probe the stacked 2-80 Hz LFP or selected band's envelope
%     (left) and the bout's stretch as a line through the top 3 PCs, the bout's
%     own or shared across all bouts (right). Yellow = the bout, light blue =
%     the other probe's own bouts. Keys: left/right = bout, v = probe 1 / 2 /
%     all, up/down = shank/depth colours, b = band, p = own/shared PCs,
%     e = LFP/envelope traces
%   vpc (cell, probes x bands) - (struct) shared PCs of all bouts used by the
%     bout viewer: coeff (channels x 3), mu (1 x channels), expl (% variance per PC)
%
% NOTE: the bout viewer's key and drawing code are local functions at the very
%   end of this script (oscpr_keypress, oscpr_redraw, oscpr_drawprobe). MATLAB
%   only sees them when the script is saved as a file, so run it from the saved file.
%
% Dependencies: NNx_BaseAnalysis.m (builds ProbeInfo),
%   NNx_editProbeInfo_ChRemove.m (Ch_Remove), oshkosh_detect_funct.m,
%   oscillation_plotEvents.m, uipickfiles, Signal Processing Toolbox (butter,
%   zp2sos, filtfilt, hilbert), Statistics Toolbox (pca, signrank, boxchart
%   needs R2020a+).
%
% NOTE (inferred data contracts, not verified on real recordings):
%   - oscillations.stim_times holds EVERY TTL-high sample, not onsets, so
%     onsets are rebuilt here from gaps in that list.
%   - prbsnippets.event_data{p}{k} rows are in oscillations.chanals order, and
%     its columns line up with oscillations.prbEvent_ext_timestamps{p}(k,:).
%     Every probe's snippet holds ALL analyzed channels (oshkosh_detect_funct.m,
%     "ALL analyzed channels"), so the bout viewer shows the other probes from
%     the bout's own snippet, at exactly the same samples. Their envelopes are
%     made when drawn, with the same filters, bins and SDs (sdch) as envz.
%   - the snippets are already bandpassed at fq_range by oshkosh (bandpass with
%     'ImpulseResponse','iir', taken as zero phase), on the whole recording, so
%     the only new edge wobble is from this script's band filter + Hilbert on
%     each snippet. edge_ms keeps that out of every window.
%   - oshkosh_detect_funct only keeps signal around events (event +/- 1s), so a
%     stim window only gets a PR if the stim onset lands inside an event AND that
%     event's snippet covers all of the windows (plus edge_ms). Stims outside an
%     event are left out, so pre/post PR (and the stim-aligned PC plot) describes
%     stims during an oscillation, not all stims. n_used tells you how many that
%     was. The sliding stim-aligned PR still takes any snippet covering each step.
%   - the stim-aligned PC plot with align_pcfit = 0 uses one PC space for every
%     recording, which assumes the channels' activity patterns don't drift much
%     across the session. If they do, align_pcfit = 1 fits each recording alone.
%   - per-window PR depends on window length AND band width: a narrow band's
%     envelope moves slowly, so a short window holds few independent points
%     (4-12 Hz in 500 ms: about 4-8). Pooled PR doesn't have that cap. Event PR
%     uses each event's own length, so check dur_ms / n_bins before comparing.
%   - 50 Hz line noise isn't notched, so it sits inside the 30-80 Hz band.
%   - oshkosh_detect_funct normalizes the TTL to its own max, so a recording
%     with no stim can read noise as stim (seen on synthetic data). 1ms pulses
%     can also be missed, since it keeps only 1 of every 30 samples.


%% NEEDS INPUT FIRST: SET UP ANALYSIS INFO

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Who is the data?
animal = '20260423-p12';

%Where is the ProbeInfo folder?
save_directory = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%Where is the raw data?
load_directory = 'E:\Roy\Silicon Probe Raw Data'; %don't change unless you're not Roy and have your data elsewhere

fq_range = [2 80]; %frequency band to check for oscillations in, Hz
probes = [1]; %which probes in the recording to look at?

%because you messed up the channel order for all pre 06-2024 animals
reorder = 0;
    % 0 = channels are fine
    % 1 = remap channels (load the matching re_index file first so newchannelidx exists)
%load 'pre202406_64_32_channel_re_index'
%load 'pre202406_64_32_channel_re_index_v2'
%load '202408_64_64_channel_re_index_v5.mat'
%load 20240708_64_32_channel_re_index.mat

%%%Analysis details%%%

%which recordings have no real stim? any condition name containing one of these words
nostim_conds = {'spon'}; %oshkosh can read TTL noise as stims here, so these get wiped below

%how long after a stim does an event still count as "with stim"?
stim_reach = 3000; %ms, 3000 matches the grey shading in the viewer. 0 = only events overlapping the TTL itself

%how close can TTL pulses be and still count as one stim (like a whisker train)?
train_gap = 1000; %ms, keep this >= post_win(2) - pre_win(1) so windows can't overlap

%what windows around stim onset do you want?
pre_win = [-500 0]; %ms
post_win = [0 500]; %ms

%what should the participation ratio be built from?
%pr_matrix = 0; %commented out: PR is now always the covariance of z-scored Hilbert envelopes (like NNx_ManifoldPCA_LvsLW.m). Correlation per window took out amplitude differences between windows
    % 0 = covariance (big-amplitude channels count for more)
    % 1 = correlation (every channel counts the same)

%%%Participation ratio details%%%
%which bands? one Hilbert amplitude envelope (and one PR) per row. Has to sit inside fq_range
hbert_bands = [4 12; 12 30; 30 80]; %Hz (theta/alpha, beta, gamma), same as NNx_ManifoldPCA_LvsLW.m
bin_ms = 5; %ms per time bin (whole number), the envelope is averaged within each bin

%z-score each channel? (without it, the loudest few channels decide the PR)
znorm = 1;
    % 0 = no
    % 1 = divide each channel by its SD over all event bins (both stim groups together, so a gain change between them isn't normalised away)

edge_ms = 300; %ms, a window is only used if its snippet covers it with this much to spare on each side (filter and Hilbert edge wobble)
n_perm = 1000; %label swaps for the pooled pre vs post test (lower it if you're in a hurry)
rng_seed = 1; %same swaps every run, so results repeat

%%%Sliding PR%%%
slide_ms = 200; %ms, length of the sliding window
slide_step_ms = 25; %ms, how far it moves each step
slide_stim_ms = [pre_win(1) post_win(2)]; %ms re: stim onset, the stretch it slides over
slide_event_ms = [-500 1000]; %ms re: event start, the stretch it slides over (snippets only reach 1 s either side of the event)

%%%PCA plot%%%
%want the events drawn as lines through PC space?
plot_pca = 1;
    % 0 = no
    % 1 = yes
pca_dims = 3;
    % 2 = PC1 vs PC2
    % 3 = PC1 vs PC2 vs PC3
ev_plot_ms = 500; %ms after event start, the stretch averaged for the event means in the bottom panel
fade_start = 0.75; %how washed out each trace starts (0 = full colour the whole way, 1 = starts white). Faded = early, full colour = late
trial_alpha = 0.3; %how see-through the single event lines are (0 = invisible, 1 = solid)

%%%Stim-aligned PCA plot%%%
%want each recording's stims drawn through PC space, lined up on stim onset?
plot_stimalign = 1;
    % 0 = no
    % 1 = yes
align_pcs = [1 2 3]; %which PCs to plot, 2 of them = flat plot, 3 = 3D (e.g. [2 3 4] to look past PC1)
align_win = [-300 300]; %ms re: stim onset
align_alpha = 0.2; %how see-through each stim's line is (0 = invisible, 1 = solid)
align_width = .7; %how thick each stim's line is
align_shift = 1;
    % 0 = leave them where they are (do they end up in the same place after the stim?)
    % 1 = move each so its onset point sits at 0 (does the stim push them the same way?)
align_pcfit = 0;
    % 0 = one PCA per probe, all recordings (same axes in every figure)
    % 1 = one PCA per recording (fits better if things drift across the session, but the axes differ between figures)

%which events should the scrollable viewer show?
view_which = 0;
    % 0 = all events
    % 1 = only with stim
    % 2 = only without stim

%want the bout viewer with each bout's PR over time and its manifold?
view_pr = 1;
    % 0 = no
    % 1 = yes (uses view_which too)

%%%nothing else to fill in%%%


%% Data Extraction: Pick recordings and detect oscillations
load(fullfile([save_directory '\' animal '-ProbeInfo.mat']))
fs = 1000; %oshkosh_detect_funct downsamples everything to 1kHz, so 1 sample = 1 ms

data_directory = uipickfiles('FilterSpec',load_directory,'Prompt', 'Choose Folders Containing Recordings you want to use [IMMEDIATE SUBFOLDERS SHOULD BE -Record Node ###-]');
data_directory = data_directory(:); %oshkosh counts files with size(x,1), so the list has to be a column or it only runs the first one

%this also makes the "all events" figure for each recording
[oscillations, snippets, prbsnippets] = oshkosh_detect_funct(animal, data_directory, probes, fq_range, save_directory, reorder, ProbeInfo);


%% Find Stim Timepoints: Rebuild stim onsets for each recording
nfiles = numel(data_directory);
file_condition = cell(nfiles,1);
n_ttl_samples = zeros(nfiles,1);
n_onsets = zeros(nfiles,1);
stim_file = [];
stim_onset = [];
for file = 1:nfiles
    %condition name, same way oshkosh gets it: everything after the last <animal>-
    occurrences = strfind(data_directory{file}, [animal '-']);
    file_condition{file} = data_directory{file}(occurrences(end)+numel(animal)+1:end);

    ttl_high = oscillations.stim_times{file}(:); %every sample where the TTL was high [samples x 1]
    n_ttl_samples(file) = numel(ttl_high);

    %no real stim in this recording, so wipe anything oshkosh found on the TTL channel
    if any(contains(file_condition{file},nostim_conds))
        ttl_high = [];
        oscillations.stim_times{file} = [];
        for prb = 1:length(probes)
            oscillations.prbEvent_during_stim{prb}(oscillations.prbEvent_file_index{prb} == file) = false;
        end
        for ch = 1:length(oscillations.chanals)
            oscillations.chEvent_during_stim{ch}(oscillations.chEvent_file_index{ch} == file) = false;
        end
    end

    if isempty(ttl_high) == false
        onsets = ttl_high([true; diff(ttl_high) > train_gap]); %first sample of each stim
        n_onsets(file) = numel(onsets);
        stim_file = [stim_file; repmat(file,numel(onsets),1)];
        stim_onset = [stim_onset; onsets];
    end
end
stim_condition = file_condition(stim_file);
stim_condition = stim_condition(:); %column, even when there are no stims
nstim = numel(stim_onset);
fileReport = table((1:nfiles)', file_condition, n_ttl_samples, n_onsets, ...
    'VariableNames', {'file','condition','n_ttl_samples','n_onsets'}); %just take a looksee after the run

%how to get PR from a [time x channels] chunk of data
%commented out: pr_matrix is gone, every PR is now the covariance of z-scored Hilbert envelopes
%if pr_matrix == 0
%    chan_matrix = @cov;
%else
%    chan_matrix = @corrcoef;
%end
%for a symmetric matrix, sum of eigenvalues = trace and sum of squared eigenvalues = sum of all squared entries, so no eig needed
pr_calc = @(C) trace(C)^2 / sum(C(:).^2);

%one bandpass per Hilbert band
nband = size(hbert_bands,1); %how many bands
signames = "LFP_hbert_" + hbert_bands(:,1)' + "to" + hbert_bands(:,2)'; %e.g. "LFP_hbert_4to12"
bandtags = hbert_bands(:,1)' + "to" + hbert_bands(:,2)'; %e.g. "4to12", for table column names
hbert_sos = cell(1,nband); %{band} second-order sections
hbert_g = zeros(1,nband); %gain of each
for j = 1:nband
    [z,p,k] = butter(2, hbert_bands(j,:)/(fs/2), 'bandpass'); %2nd order per edge, same as NNx_ManifoldPCA_LvsLW.m
    [hbert_sos{j},hbert_g(j)] = zp2sos(z,p,k); %second-order sections, stable at low cutoffs
end

stim_conds = unique(stim_condition,'stable'); %stim conditions, in the order they came in
stim_conds = stim_conds(:);
nstimcond = numel(stim_conds);
group_names = {'no stim'; 'with stim'};
rng(rng_seed);
run_notes = strings(0,1); %problems for the report at the end


%% Data PreProcessing: Hilbert envelope of every event snippet, binned and z-scored
envz = cell(1,length(probes)); %{probe}{band}{event} [bins x channels]
binstart = cell(1,length(probes)); %{probe}{event} [bins x 1] recording sample where each bin starts
evmask = cell(1,length(probes)); %{probe}{event} [bins x 1] true = inside the event itself
sdch = cell(1,length(probes)); %{probe} [bands x channels] SD used for the z-score
nch_prb = zeros(1,length(probes)); %channels on each probe

tic
for prb = 1:length(probes)
    ev_times = oscillations.prbEvent_timestamps{prb}; %[events x 2], start and stop samples
    ev_ext = oscillations.prbEvent_ext_timestamps{prb}; %[events x 2], same but +/- 1s, what the snippets cover
    prb_chans = find(oscillations.chan_probe == prb); %snippets hold every analyzed channel, so this picks out this probe's rows
    nev = size(ev_times,1);
    nch_prb(prb) = numel(prb_chans);
    disp(['PR: Hilbert envelopes for probe ' num2str(probes(prb)) ', ' num2str(nev) ' events'])

    envz{prb} = cell(1,nband);
    for j = 1:nband
        envz{prb}{j} = cell(nev,1);
    end
    binstart{prb} = cell(nev,1);
    evmask{prb} = cell(nev,1);
    for ev = 1:nev
        snip = prbsnippets.event_data{prb}{ev}(prb_chans,:)'; %[snippet samples x channels], already fq_range bandpassed by oshkosh
        nbin = floor(size(snip,1)/bin_ms);
        binstart{prb}{ev} = ev_ext(ev,1) + bin_ms*(0:nbin-1)'; %recording sample where each bin starts
        evmask{prb}{ev} = binstart{prb}{ev} >= ev_times(ev,1) & binstart{prb}{ev} + bin_ms - 1 <= ev_times(ev,2); %bins inside the event
        for j = 1:nband
            env = abs(hilbert(filtfilt(hbert_sos{j}, hbert_g(j), snip))); %[snippet samples x channels] amplitude envelope of this band, uV
            env = env(1:nbin*bin_ms,:); %drop the leftover samples that don't fill a bin
            envz{prb}{j}{ev} = single(reshape(mean(reshape(env,bin_ms,nbin,[]),1),nbin,[])); %[bins x channels] mean per bin, single halves the memory
        end
    end

    %SD of each channel over the bins inside every event (both stim groups, all files together)
    sdch{prb} = ones(nband,nch_prb(prb));
    if znorm == 1
        for j = 1:nband
            evpool = [];
            for ev = 1:nev
                evpool = [evpool; envz{prb}{j}{ev}(evmask{prb}{ev},:)]; %#ok<AGROW>
            end
            sdch{prb}(j,:) = std(double(evpool),0,1);
            for ev = 1:nev
                envz{prb}{j}{ev} = envz{prb}{j}{ev} ./ single(sdch{prb}(j,:)); %z-units
            end
        end
        clear evpool
    end
end
clear snip env
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data PreProcessing: Sort probe events by stim present or not, and get each event's PR
eventTable = cell(1,length(probes));
stimSortSummary = cell(1,length(probes));
with_stim = cell(1,length(probes));
for prb = 1:length(probes)
    ev_times = oscillations.prbEvent_timestamps{prb}; %[events x 2], start and stop samples
    %ev_ext = oscillations.prbEvent_ext_timestamps{prb}; %only the old raw-sample PR needed this
    ev_file = oscillations.prbEvent_file_index{prb};
    during_stim = oscillations.prbEvent_during_stim{prb};
    %prb_chans = find(oscillations.chan_probe == prb); %only the old raw-sample PR needed this
    nev = size(ev_times,1);
    nch = nch_prb(prb);

    with_stim{prb} = false(nev,1);
    %pr = zeros(nev,1); %old: one PR per event on the raw 1 kHz LFP
    pr = nan(nev,nband); %per-window PR of each event, one column per band
    n_bins = zeros(nev,1);
    for ev = 1:nev
        %does a stim overlap this event, or come just before it?
        ttl_high = oscillations.stim_times{ev_file(ev)};
        with_stim{prb}(ev) = during_stim(ev) == 1 | any(ttl_high <= ev_times(ev,1) & ttl_high >= ev_times(ev,1)-stim_reach);

        %PR over the event itself, without the extra 1s on each side
        %old: raw 1 kHz fq_range LFP, covariance or correlation. Swapped for the binned z-scored Hilbert envelope of each band
        %ev_idx = (ev_times(ev,1):ev_times(ev,2)) - ev_ext(ev,1) + 1; %recording samples to snippet samples
        %ev_dat = prbsnippets.event_data{prb}{ev}(prb_chans,ev_idx); %[channels x time]
        %pr(ev) = pr_calc(chan_matrix(ev_dat'));
        n_bins(ev) = nnz(evmask{prb}{ev});
        for j = 1:nband
            pr(ev,j) = pr_calc(cov(double(envz{prb}{j}{ev}(evmask{prb}{ev},:)))); %[event bins x channels]
        end
    end

    eventTable{prb} = table(ev_file, oscillations.prbEvent_condition{prb}, ev_times(:,1), ev_times(:,2), ...
        oscillations.prbEvent_durations{prb}, oscillations.prbEvent_peakNormedPower{prb}, ...
        logical(during_stim), with_stim{prb}, n_bins, ...
        'VariableNames', {'file','condition','start','stop','dur_ms','peak_rms','during_stim','with_stim','n_bins'});
    for j = 1:nband
        eventTable{prb}.(char("pr_" + bandtags(j))) = pr(:,j);
    end

    %medians for each group, and the pooled PR of each group
    n_events = zeros(2,1);
    med_dur_ms = zeros(2,1);
    med_peak_rms = zeros(2,1);
    %med_pr = zeros(2,1); %old: one PR column, now one per band
    med_pr = nan(2,nband);
    pool_pr = nan(2,nband);
    pooleq_pr = nan(2,nband);
    for g = 1:2
        in_group = eventTable{prb}.with_stim == (g == 2);
        n_events(g) = sum(in_group);
        med_dur_ms(g) = median(eventTable{prb}.dur_ms(in_group));
        med_peak_rms(g) = median(eventTable{prb}.peak_rms(in_group));
        for j = 1:nband
            med_pr(g,j) = median(pr(in_group,j));
            %pooled: each event mean-centred, then all of the group's event bins in one covariance
            %equal-weight pooled: same, but each event scaled to the same total variance first, so loud events don't dominate
            S = zeros(nch,nch);
            Seq = zeros(nch,nch);
            for ev = find(in_group)'
                X = double(envz{prb}{j}{ev}(evmask{prb}{ev},:)); %[event bins x channels]
                X = X - mean(X,1);
                XX = X'*X;
                S = S + XX; %PR doesn't care about the scale, so the plain sum is fine
                Seq = Seq + XX/trace(XX); %trace = this event's total variance
            end
            pool_pr(g,j) = pr_calc(S); %NaN when the group has no events
            pooleq_pr(g,j) = pr_calc(Seq);
        end
    end
    stimSortSummary{prb} = table(group_names, n_events, med_dur_ms, med_peak_rms, ...
        'VariableNames', {'group','n_events','med_dur_ms','med_peak_rms'});
    for j = 1:nband
        stimSortSummary{prb}.(char("med_pr_" + bandtags(j))) = med_pr(:,j);
        stimSortSummary{prb}.(char("poolpr_" + bandtags(j))) = pool_pr(:,j);
        stimSortSummary{prb}.(char("pooleqpr_" + bandtags(j))) = pooleq_pr(:,j);
    end
end


%% Data PreProcessing: Participation ratio in pre and post stim windows (per-window and pooled)
%commented out: old version, raw 1 kHz LFP, one PR per window, any event covering each window on its own
%{
pr_pre = cell(1,length(probes));
pr_post = cell(1,length(probes));
prepostSummary = cell(1,length(probes));
stim_conds = unique(stim_condition,'stable');
for prb = 1:length(probes)
    ev_ext = oscillations.prbEvent_ext_timestamps{prb};
    ev_file = oscillations.prbEvent_file_index{prb};
    prb_chans = find(oscillations.chan_probe == prb);
    pr_pre{prb} = NaN(nstim,1); %[stims x 1], stays NaN if no snippet covers the window
    pr_post{prb} = NaN(nstim,1);

    for s = 1:nstim
        for w = 1:2 %1 = pre, 2 = post
            if w == 1
                win = stim_onset(s) + pre_win;
            else
                win = stim_onset(s) + post_win;
            end
            win_samples = win(1):win(2)-1; %win(2)-1 so pre and post don't share the onset sample

            %timestamps restart in each recording, so only look in this one
            ev = find(ev_file == stim_file(s) & ev_ext(:,1) <= win_samples(1) & ev_ext(:,2) >= win_samples(end),1);
            if isempty(ev) == false
                win_dat = prbsnippets.event_data{prb}{ev}(prb_chans,win_samples-ev_ext(ev,1)+1); %[channels x time]
                if w == 1
                    pr_pre{prb}(s) = pr_calc(chan_matrix(win_dat'));
                else
                    pr_post{prb}(s) = pr_calc(chan_matrix(win_dat'));
                end
            end
        end
    end

    %split by condition so whisker and light don't get mushed together
    ncond = numel(stim_conds);
    n_stim = zeros(ncond,1);
    n_used = zeros(ncond,1);
    mean_pr_pre = zeros(ncond,1);
    mean_pr_post = zeros(ncond,1);
    for cnd = 1:ncond
        in_cond = strcmp(stim_condition,stim_conds{cnd});
        both_win = in_cond & ~isnan(pr_pre{prb}) & ~isnan(pr_post{prb}); %only stims with both windows, so pre vs post is paired
        n_stim(cnd) = sum(in_cond);
        n_used(cnd) = sum(both_win);
        mean_pr_pre(cnd) = mean(pr_pre{prb}(both_win));
        mean_pr_post(cnd) = mean(pr_post{prb}(both_win));
    end
    prepostSummary{prb} = table(stim_conds, n_stim, n_used, mean_pr_pre, mean_pr_post, ...
        'VariableNames', {'condition','n_stim','n_used','pr_pre','pr_post'});
end
%}

%per-window: one PR per stim window, then the mean over stims. Capped by how few independent points a window holds
%pooled: every stim's window mean-centred, all of them in one covariance, one PR. Tested by swapping pre and post within random stims
pr_pre = cell(1,length(probes));
pr_post = cell(1,length(probes));
stim_ev = cell(1,length(probes));
prepostSummary = cell(1,length(probes));
winlen_s = (pre_win(2)-pre_win(1))/1000; %s, length of one window (pre and post should match)

tic
for prb = 1:length(probes)
    ev_times = oscillations.prbEvent_timestamps{prb}; %[events x 2], start and stop samples
    ev_ext = oscillations.prbEvent_ext_timestamps{prb};
    ev_file = oscillations.prbEvent_file_index{prb};
    nch = nch_prb(prb);
    pr_pre{prb} = nan(nstim,nband); %[stims x bands], stays NaN if no snippet covers the windows
    pr_post{prb} = nan(nstim,nband);
    stim_ev{prb} = nan(nstim,1);
    Spre = cell(nstim,nband); %{stim x band} [channels x channels] mean-centred scatter of each window, for the pooled PR
    Spost = cell(nstim,nband);

    for s = 1:nstim
        %one event's snippet has to cover pre AND post (plus edge_ms), so the pair is from the same stretch of data
        a = stim_onset(s) + pre_win(1);
        b = stim_onset(s) + post_win(2);
        %ev = find(ev_file == stim_file(s) & ev_ext(:,1) <= a - edge_ms & ev_ext(:,2) >= b - 1 + edge_ms, 1); %old: any snippet (event +/- 1 s) covering the windows counted, even if the stim missed the event itself
        ev = find(ev_file == stim_file(s) & ev_times(:,1) <= stim_onset(s) & ev_times(:,2) >= stim_onset(s) ...
            & ev_ext(:,1) <= a - edge_ms & ev_ext(:,2) >= b - 1 + edge_ms, 1); %stim onset inside the event, and its snippet covers pre + post + edge_ms. Timestamps restart in each recording, so only look in this one
        if isempty(ev) == true
            continue
        end
        stim_ev{prb}(s) = ev;
        sb = binstart{prb}{ev};
        inpre = sb >= stim_onset(s) + pre_win(1) & sb + bin_ms - 1 <= stim_onset(s) + pre_win(2) - 1; %pre bins end before the onset sample
        inpost = sb >= stim_onset(s) + post_win(1) & sb + bin_ms - 1 <= stim_onset(s) + post_win(2) - 1;
        for j = 1:nband
            Xa = double(envz{prb}{j}{ev}(inpre,:)); %[pre bins x channels]
            Xb = double(envz{prb}{j}{ev}(inpost,:)); %[post bins x channels]
            pr_pre{prb}(s,j) = pr_calc(cov(Xa));
            pr_post{prb}(s,j) = pr_calc(cov(Xb));
            Xa = Xa - mean(Xa,1);
            Xb = Xb - mean(Xb,1);
            Spre{s,j} = Xa'*Xa;
            Spost{s,j} = Xb'*Xb;
        end
    end

    %split by condition so whisker and light don't get mushed together
    n_stim = zeros(nstimcond,1);
    n_used = zeros(nstimcond,1);
    prwin_pre = nan(nstimcond,nband);
    prwin_post = nan(nstimcond,nband);
    sr_p = nan(nstimcond,nband);
    prpool_pre = nan(nstimcond,nband);
    prpool_post = nan(nstimcond,nband);
    perm_p = nan(nstimcond,nband);
    prpooleq_pre = nan(nstimcond,nband);
    prpooleq_post = nan(nstimcond,nband);
    permeq_p = nan(nstimcond,nband);
    n_indep = nan(nstimcond,nband);
    for cnd = 1:nstimcond
        in_cond = strcmp(stim_condition,stim_conds{cnd});
        idx = find(in_cond & isnan(stim_ev{prb}) == false); %stims with both windows, so pre vs post is paired
        n_stim(cnd) = sum(in_cond);
        n_used(cnd) = numel(idx);
        if isempty(idx) == true
            run_notes(end+1) = "Probe " + probes(prb) + ", " + stim_conds{cnd} + ": no stim landed inside an event whose snippet covers pre + post + edge_ms"; %#ok<SAGROW>
            continue
        end
        for j = 1:nband
            prwin_pre(cnd,j) = mean(pr_pre{prb}(idx,j));
            prwin_post(cnd,j) = mean(pr_post{prb}(idx,j));
            sr_p(cnd,j) = signrank(pr_post{prb}(idx,j),pr_pre{prb}(idx,j));

            SA = cat(3,Spre{idx,j}); %[channels x channels x stims]
            SB = cat(3,Spost{idx,j});
            prpool_pre(cnd,j) = pr_calc(sum(SA,3));
            prpool_post(cnd,j) = pr_calc(sum(SB,3));
            realdiff = prpool_post(cnd,j) - prpool_pre(cnd,j);

            %swap pre and post within random stims, so the pairing stays
            nulldiff = zeros(n_perm,1);
            for it = 1:n_perm
                swap = rand(1,1,numel(idx)) < 0.5; %[1 x 1 x stims] true = swap this stim's pre and post
                A = sum(SA.*(swap == false) + SB.*swap,3);
                B = sum(SB.*(swap == false) + SA.*swap,3);
                nulldiff(it) = pr_calc(B) - pr_calc(A);
            end
            perm_p(cnd,j) = (sum(abs(nulldiff) >= abs(realdiff))+1)/(n_perm+1); %two-sided

            %equal-weight pooled: each window scaled to the same total variance (trace = 1) before stacking, so loud windows don't dominate
            trA = sum(SA.*eye(nch),[1 2]); %[1 x 1 x stims] total variance of each pre window
            trB = sum(SB.*eye(nch),[1 2]);
            SAn = SA./trA;
            SBn = SB./trB;
            prpooleq_pre(cnd,j) = pr_calc(sum(SAn,3));
            prpooleq_post(cnd,j) = pr_calc(sum(SBn,3));
            realdiff = prpooleq_post(cnd,j) - prpooleq_pre(cnd,j);
            nulldiff = zeros(n_perm,1);
            for it = 1:n_perm
                swap = rand(1,1,numel(idx)) < 0.5;
                A = sum(SAn.*(swap == false) + SBn.*swap,3);
                B = sum(SBn.*(swap == false) + SAn.*swap,3);
                nulldiff(it) = pr_calc(B) - pr_calc(A);
            end
            permeq_p(cnd,j) = (sum(abs(nulldiff) >= abs(realdiff))+1)/(n_perm+1); %two-sided
            n_indep(cnd,j) = n_used(cnd)*winlen_s*(hbert_bands(j,2)-hbert_bands(j,1)); %rough: an envelope of a band B Hz wide changes about B times a second
        end
    end
    prepostSummary{prb} = table(stim_conds, n_stim, n_used, 'VariableNames', {'condition','n_stim','n_used'});
    for j = 1:nband
        prepostSummary{prb}.(char("prwin_pre_" + bandtags(j))) = prwin_pre(:,j);
        prepostSummary{prb}.(char("prwin_post_" + bandtags(j))) = prwin_post(:,j);
        prepostSummary{prb}.(char("signrank_p_" + bandtags(j))) = sr_p(:,j);
        prepostSummary{prb}.(char("prpool_pre_" + bandtags(j))) = prpool_pre(:,j);
        prepostSummary{prb}.(char("prpool_post_" + bandtags(j))) = prpool_post(:,j);
        prepostSummary{prb}.(char("perm_p_" + bandtags(j))) = perm_p(:,j);
        prepostSummary{prb}.(char("prpooleq_pre_" + bandtags(j))) = prpooleq_pre(:,j);
        prepostSummary{prb}.(char("prpooleq_post_" + bandtags(j))) = prpooleq_post(:,j);
        prepostSummary{prb}.(char("permeq_p_" + bandtags(j))) = permeq_p(:,j);
        prepostSummary{prb}.(char("n_indep_" + bandtags(j))) = n_indep(:,j);
    end
end
clear Spre Spost SA SB SAn SBn A B
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data PreProcessing: Sliding window PR over time, per-window and pooled
%stim-aligned (per stim condition) and event-aligned (no stim / with stim). Every step uses the same window
%length, so the steps can be compared with each other. Pooled = all the stims/events at that step in one covariance
stepst_stim = slide_stim_ms(1):slide_step_ms:(slide_stim_ms(2)-slide_ms); %window starts, ms re: stim onset
slide_t_stim = stepst_stim + slide_ms/2; %window centres
stepst_ev = slide_event_ms(1):slide_step_ms:(slide_event_ms(2)-slide_ms); %window starts, ms re: event start
slide_t_ev = stepst_ev + slide_ms/2;
nstep_stim = numel(stepst_stim);
nstep_ev = numel(stepst_ev);

slidePR_stim_each = cell(1,length(probes));
slidePR_stim_pool = cell(1,length(probes));
slidePR_stim_pooleq = cell(1,length(probes));
slide_n_stim = cell(1,length(probes));
slidePR_ev_each = cell(1,length(probes));
slidePR_ev_pool = cell(1,length(probes));
slidePR_ev_pooleq = cell(1,length(probes));
slide_n_ev = cell(1,length(probes));

tic
for prb = 1:length(probes)
    disp(['PR: Sliding windows for probe ' num2str(probes(prb))])
    ev_times = oscillations.prbEvent_timestamps{prb};
    ev_ext = oscillations.prbEvent_ext_timestamps{prb};
    ev_file = oscillations.prbEvent_file_index{prb};
    nev = size(ev_times,1);
    nch = nch_prb(prb);

    %stim-aligned
    slidePR_stim_each{prb} = nan(nstim,nstep_stim,nband);
    slidePR_stim_pool{prb} = nan(nstimcond,nstep_stim,nband);
    slidePR_stim_pooleq{prb} = nan(nstimcond,nstep_stim,nband);
    slide_n_stim{prb} = zeros(nstimcond,nstep_stim);
    for cnd = 1:nstimcond
        for t = 1:nstep_stim
            S = zeros(nch,nch,nband);
            Seq = zeros(nch,nch,nband); %equal-weight: each window scaled to the same total variance
            for s = find(strcmp(stim_condition,stim_conds{cnd}))'
                a = stim_onset(s) + stepst_stim(t);
                b = a + slide_ms;
                ev = find(ev_file == stim_file(s) & ev_ext(:,1) <= a - edge_ms & ev_ext(:,2) >= b - 1 + edge_ms, 1); %an event snippet that covers this step
                if isempty(ev) == true
                    continue
                end
                sb = binstart{prb}{ev};
                inw = sb >= a & sb + bin_ms - 1 <= b - 1; %bins in this window
                slide_n_stim{prb}(cnd,t) = slide_n_stim{prb}(cnd,t) + 1;
                for j = 1:nband
                    X = double(envz{prb}{j}{ev}(inw,:)); %[window bins x channels]
                    slidePR_stim_each{prb}(s,t,j) = pr_calc(cov(X));
                    X = X - mean(X,1);
                    XX = X'*X;
                    S(:,:,j) = S(:,:,j) + XX;
                    Seq(:,:,j) = Seq(:,:,j) + XX/trace(XX);
                end
            end
            for j = 1:nband
                slidePR_stim_pool{prb}(cnd,t,j) = pr_calc(S(:,:,j)); %NaN when no stim covers this step
                slidePR_stim_pooleq{prb}(cnd,t,j) = pr_calc(Seq(:,:,j));
            end
        end
    end

    %event-aligned, each event's own snippet (it reaches 1 s either side of the event)
    slidePR_ev_each{prb} = nan(nev,nstep_ev,nband);
    slidePR_ev_pool{prb} = nan(2,nstep_ev,nband);
    slidePR_ev_pooleq{prb} = nan(2,nstep_ev,nband);
    slide_n_ev{prb} = zeros(2,nstep_ev);
    for g = 1:2
        evs = find(with_stim{prb} == (g == 2))';
        for t = 1:nstep_ev
            S = zeros(nch,nch,nband);
            Seq = zeros(nch,nch,nband);
            for ev = evs
                a = ev_times(ev,1) + stepst_ev(t);
                b = a + slide_ms;
                if ev_ext(ev,1) > a - edge_ms || ev_ext(ev,2) < b - 1 + edge_ms
                    continue %the snippet doesn't reach this step
                end
                sb = binstart{prb}{ev};
                inw = sb >= a & sb + bin_ms - 1 <= b - 1;
                slide_n_ev{prb}(g,t) = slide_n_ev{prb}(g,t) + 1;
                for j = 1:nband
                    X = double(envz{prb}{j}{ev}(inw,:));
                    slidePR_ev_each{prb}(ev,t,j) = pr_calc(cov(X));
                    X = X - mean(X,1);
                    XX = X'*X;
                    S(:,:,j) = S(:,:,j) + XX;
                    Seq(:,:,j) = Seq(:,:,j) + XX/trace(XX);
                end
            end
            for j = 1:nband
                slidePR_ev_pool{prb}(g,t,j) = pr_calc(S(:,:,j));
                slidePR_ev_pooleq{prb}(g,t,j) = pr_calc(Seq(:,:,j));
            end
        end
    end
end
clear S Seq X XX
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);


%% Data PreProcessing: How many events go with each stim?
%you'd expect one event per stim. An event is linked to a stim if it starts within stim_reach after the onset, or spans it
stimEventReport = cell(1,length(probes));
for prb = 1:length(probes)
    ev_times = oscillations.prbEvent_timestamps{prb};
    ev_file = oscillations.prbEvent_file_index{prb};
    stimEventReport{prb} = table();
    for file = 1:nfiles
        ons = stim_onset(stim_file == file);
        if isempty(ons) == true
            continue %no stims in this recording
        end
        evf = find(ev_file == file);
        link = false(numel(ons),numel(evf)); %[stims x events] true = linked
        for s = 1:numel(ons)
            link(s,:) = (ev_times(evf,1) >= ons(s) & ev_times(evf,1) <= ons(s) + stim_reach)' | ...
                (ev_times(evf,1) <= ons(s) & ev_times(evf,2) >= ons(s))';
        end
        nlink = sum(link,2); %events per stim
        stimEventReport{prb} = [stimEventReport{prb}; table(file, string(file_condition{file}), numel(ons), sum(nlink == 0), sum(nlink == 1), sum(nlink >= 2), ...
            sum(any(link,1)), sum(sum(link,1) >= 2), ...
            'VariableNames', {'file','condition','n_stims','stims_0ev','stims_1ev','stims_2plusev','events_linked','events_multistim'})];
    end
end


%% Data Visualization: Events by stim group and pre vs post participation ratio
%commented out: old single-band figure (one PR column). Replaced by the per-band figures below
%{
for prb = 1:length(probes)
    figure('Name',[animal ' probe ' num2str(probes(prb)) ' oscillation events'],'Color','w');
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(1,3,'TileSpacing','compact');
    title(tl,[animal ' - probe ' num2str(probes(prb))],'Interpreter','none');

    %every probe event, split by stim group
    stim_group = repmat({'no stim'},height(eventTable{prb}),1);
    stim_group(eventTable{prb}.with_stim) = {'with stim'};
    stim_group = categorical(stim_group,{'no stim','with stim'});

    ax = nexttile(tl);
    boxchart(ax,stim_group,eventTable{prb}.dur_ms);
    ylabel(ax,'Event duration (ms)');

    ax = nexttile(tl);
    boxchart(ax,stim_group,eventTable{prb}.pr);
    ylabel(ax,'Participation ratio over event');
    ylim(ax,[1 inf]); %1 is the lowest PR can go

    %pre vs post stim, one line per condition. mean +/- SEM across stims with both windows
    ax = nexttile(tl);
    hold(ax,'on');
    for cnd = 1:numel(stim_conds)
        both_win = strcmp(stim_condition,stim_conds{cnd}) & ~isnan(pr_pre{prb}) & ~isnan(pr_post{prb});
        means = [mean(pr_pre{prb}(both_win)) mean(pr_post{prb}(both_win))];
        sems = [std(pr_pre{prb}(both_win)) std(pr_post{prb}(both_win))] / sqrt(sum(both_win));
        errorbar(ax,[1 2],means,sems,'-o','LineWidth',2,'MarkerSize',8,'DisplayName',[stim_conds{cnd} ' (n=' num2str(sum(both_win)) ')']);
    end
    hold(ax,'off');
    xlim(ax,[0.5 2.5]);
    ylim(ax,[1 inf]);
    xticks(ax,[1 2]);
    xticklabels(ax,{'pre','post'});
    xlabel(ax,sprintf('pre = [%d %d] ms, post = [%d %d] ms from stim onset',pre_win,post_win));
    ylabel(ax,'Participation ratio');
    legend(ax,'Interpreter','none','Location','best');
    box(ax,'off');
end
%}

grpcol = [0.5 0.5 0.5; 0.85 0.15 0.15]; %no stim = grey, with stim = red
condcol = lines(max(nstimcond,1)); %one colour per stim condition

for prb = 1:length(probes)
    %events: duration, then per-event PR (boxes) and pooled PR (diamond) for each band
    figure('Name',[animal '-probe' num2str(probes(prb)) '-events'],'Color','w');
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(1,nband+1,'TileSpacing','compact');
    title(tl,{[animal ' - probe ' num2str(probes(prb)) ': oscillation events'],'boxes = per-event PR, diamond = pooled PR of the group'},'Interpreter','none');

    stim_group = repmat({'no stim'},height(eventTable{prb}),1);
    stim_group(eventTable{prb}.with_stim) = {'with stim'};
    stim_group = categorical(stim_group,{'no stim','with stim'});

    ax = nexttile(tl);
    boxchart(ax,stim_group,eventTable{prb}.dur_ms);
    ylabel(ax,'Event duration (ms)');

    for j = 1:nband
        ax = nexttile(tl);
        boxchart(ax,stim_group,eventTable{prb}.(char("pr_" + bandtags(j))));
        hold(ax,'on');
        plot(ax,categorical(group_names',{'no stim','with stim'}),stimSortSummary{prb}.(char("poolpr_" + bandtags(j)))','kd','MarkerFaceColor','k','MarkerSize',8);
        hold(ax,'off');
        ylabel(ax,'Participation ratio over event');
        ylim(ax,[1 inf]); %1 is the lowest PR can go
        title(ax,signames(j),'Interpreter','none');
    end

    %pre vs post: top row per-window (mean +- SEM over stims), bottom row pooled
    if nstimcond > 0
        figure('Name',[animal '-probe' num2str(probes(prb)) '-prepost'],'Color','w');
        set(gcf, 'Position', get(0, 'Screensize'));
        tl = tiledlayout(3,nband,'TileSpacing','compact');
        title(tl,{[animal ' - probe ' num2str(probes(prb)) ': PR pre vs post stim'], ...
            sprintf('pre = [%d %d] ms, post = [%d %d] ms. top = per-window (mean +- SEM, signrank p), middle = pooled, bottom = equal-weight pooled (swap test p)',pre_win,post_win)},'Interpreter','none');
        for row = 1:3
            for j = 1:nband
                ax = nexttile(tl);
                hold(ax,'on');
                for cnd = 1:nstimcond
                    r = prepostSummary{prb}(cnd,:);
                    if r.n_used == 0
                        continue
                    end
                    idx = find(strcmp(stim_condition,stim_conds{cnd}) & isnan(stim_ev{prb}) == false);
                    if row == 1
                        means = [r.(char("prwin_pre_" + bandtags(j))) r.(char("prwin_post_" + bandtags(j)))];
                        sems = [std(pr_pre{prb}(idx,j)) std(pr_post{prb}(idx,j))]/sqrt(numel(idx));
                        errorbar(ax,[1 2],means,sems,'-o','Color',condcol(cnd,:),'LineWidth',2,'MarkerSize',8, ...
                            'DisplayName',[stim_conds{cnd} ' (n=' num2str(r.n_used) ', p=' num2str(r.(char("signrank_p_" + bandtags(j))),'%.3f') ')']);
                    elseif row == 2
                        plot(ax,[1 2],[r.(char("prpool_pre_" + bandtags(j))) r.(char("prpool_post_" + bandtags(j)))],'-o','Color',condcol(cnd,:),'LineWidth',2,'MarkerSize',8, ...
                            'MarkerFaceColor',condcol(cnd,:),'DisplayName',[stim_conds{cnd} ' (n=' num2str(r.n_used) ', p=' num2str(r.(char("perm_p_" + bandtags(j))),'%.3f') ')']);
                    else
                        plot(ax,[1 2],[r.(char("prpooleq_pre_" + bandtags(j))) r.(char("prpooleq_post_" + bandtags(j)))],'-s','Color',condcol(cnd,:),'LineWidth',2,'MarkerSize',8, ...
                            'MarkerFaceColor',condcol(cnd,:),'DisplayName',[stim_conds{cnd} ' (n=' num2str(r.n_used) ', p=' num2str(r.(char("permeq_p_" + bandtags(j))),'%.3f') ')']);
                    end
                end
                hold(ax,'off');
                xlim(ax,[0.5 2.5]);
                ylim(ax,[1 inf]);
                xticks(ax,[1 2]);
                xticklabels(ax,{'pre','post'});
                if row == 1
                    ylabel(ax,'PR, per-window mean');
                elseif row == 2
                    ylabel(ax,'PR, pooled');
                else
                    ylabel(ax,'PR, equal-weight pooled');
                end
                title(ax,signames(j),'Interpreter','none');
                legend(ax,'Interpreter','none','Location','best');
                box(ax,'off');
            end
        end
    end
end


%% Data Visualization: Sliding window PR over time
%rows = bands. Left = stim-aligned (one colour per condition), right = event-aligned (grey = no stim, red = with stim)
%solid + shading = per-window mean +- SEM, dashed = pooled. Thin grey steps (right axis) = how many stims/events went into each step
for prb = 1:length(probes)
    figure('Name',[animal '-probe' num2str(probes(prb)) '-slidingPR'],'Color','w');
    set(gcf, 'Position', get(0, 'Screensize'));
    tl = tiledlayout(nband,2,'TileSpacing','compact');
    title(tl,{[animal ' - probe ' num2str(probes(prb)) ': sliding PR, ' num2str(slide_ms) ' ms window, ' num2str(slide_step_ms) ' ms steps'], ...
        'solid + shading = per-window mean +- SEM, dashed = pooled, dotted = equal-weight pooled, grey steps (right axis) = n'},'Interpreter','none');
    for j = 1:nband
        for col = 1:2
            ax = nexttile(tl);
            hold(ax,'on');
            if col == 1
                tt = slide_t_stim;
                ngrp = nstimcond;
            else
                tt = slide_t_ev;
                ngrp = 2;
            end
            nall = zeros(1,numel(tt));
            for g = 1:ngrp
                if col == 1
                    rows = strcmp(stim_condition,stim_conds{g});
                    each = slidePR_stim_each{prb}(rows,:,j); %[stims x steps]
                    pool = slidePR_stim_pool{prb}(g,:,j);
                    pooleq = slidePR_stim_pooleq{prb}(g,:,j);
                    nn = slide_n_stim{prb}(g,:);
                    gc = condcol(g,:);
                    gname = stim_conds{g};
                else
                    rows = with_stim{prb} == (g == 2);
                    each = slidePR_ev_each{prb}(rows,:,j); %[events x steps]
                    pool = slidePR_ev_pool{prb}(g,:,j);
                    pooleq = slidePR_ev_pooleq{prb}(g,:,j);
                    nn = slide_n_ev{prb}(g,:);
                    gc = grpcol(g,:);
                    gname = group_names{g};
                end
                if all(nn == 0)
                    continue
                end
                mn = mean(each,1,'omitnan');
                se = std(each,0,1,'omitnan')./sqrt(nn);
                ok = isnan(mn) == false & isnan(se) == false;
                fill(ax,[tt(ok) fliplr(tt(ok))],[mn(ok)+se(ok) fliplr(mn(ok)-se(ok))],gc,'FaceAlpha',0.2,'EdgeColor','none','HandleVisibility','off');
                plot(ax,tt,mn,'-','Color',gc,'LineWidth',1.5,'DisplayName',[gname ' per-window']);
                plot(ax,tt,pool,'--','Color',gc,'LineWidth',1.5,'DisplayName',[gname ' pooled']);
                plot(ax,tt,pooleq,':','Color',gc,'LineWidth',2,'DisplayName',[gname ' equal-weight pooled']);
                nall = nall + nn;
            end
            xline(ax,0,'HandleVisibility','off');
            ylabel(ax,'PR');
            yyaxis(ax,'right');
            stairs(ax,tt,nall,'Color',[0.6 0.6 0.6],'HandleVisibility','off');
            ylabel(ax,'n');
            ax.YAxis(2).Color = [0.6 0.6 0.6];
            yyaxis(ax,'left');
            hold(ax,'off');
            if col == 1
                xlabel(ax,'ms re: stim onset (window centre)');
                title(ax,[char(signames(j)) ', stim-aligned'],'Interpreter','none');
            else
                xlabel(ax,'ms re: event start (window centre)');
                title(ax,[char(signames(j)) ', event-aligned'],'Interpreter','none');
            end
            legend(ax,'Interpreter','none','Location','best','FontSize',7);
            box(ax,'off');
        end
    end
end


%% Data Visualization: Events as lines through PC space
%one figure per probe, one column per band. One PCA per band, fit to the bins inside every event (both groups).
%top = every event start to stop, grey = no stim, red = with stim. bottom = means: each group over the first
%ev_plot_ms after event start, and each stim condition over pre_win(1) to post_win(2) around the stim (if there are stims).
%direction: faded = early, full colour = late, dot = start (event start, or stim onset), square = end
if plot_pca == 1
    disp('PCA: Making some pictures for you')
    nbp = round(ev_plot_ms/bin_ms); %bins in the event-start-aligned means
    nbs = round((post_win(2)-pre_win(1))/bin_ms); %bins in the stim-aligned means
    onset_bs = round(-pre_win(1)/bin_ms) + 1; %bin of the stim-aligned means where the stim comes on
    for prb = 1:length(probes)
        ev_times = oscillations.prbEvent_timestamps{prb};
        nev = size(ev_times,1);
        nch = nch_prb(prb);
        if nev < 2
            run_notes(end+1) = "Probe " + probes(prb) + ": no PCA plot, fewer than 2 events"; %#ok<SAGROW>
            continue
        end

        figure('Name',[animal '-probe' num2str(probes(prb)) '-eventPCA'],'Color','w');
        set(gcf, 'Position', get(0, 'Screensize'));
        tiledlayout(2,nband) %top row = all events, bottom row = means, one column per band

        for j = 1:nband
            %one PCA on the bins inside every event
            Xall = [];
            for ev = 1:nev
                Xall = [Xall; double(envz{prb}{j}{ev}(evmask{prb}{ev},:))]; %#ok<AGROW>
            end
            [coeff,~,latent,~,~,mu] = pca(Xall,'NumComponents',pca_dims); %coeff is [channels x pca_dims], mu is the mean that was taken out
            pcexpl = 100*latent/sum(latent); %percent variance per PC
            clear Xall

            %top panel: every event, thin and see-through, faded at the start. Random drawing order so no colour always sits on top
            nexttile(j)
            hold on
            for ev = randperm(nev)
                P = (double(envz{prb}{j}{ev}(evmask{prb}{ev},:)) - mu)*coeff; %[event bins x pca_dims]
                gc = grpcol(with_stim{prb}(ev)+1,:);
                fadecol = gc + fade_start*(1 - gc); %washed out version of the group colour, the start of the trace
                tracecol = fadecol + linspace(0,1,size(P,1))'*(gc - fadecol); %[event bins x RGB] faded to full colour
                %a patch with no face and interpolated edge colours is a line whose colour changes along it. The NaN stops it closing into a loop
                if pca_dims == 3
                    patch('XData',[P(:,1); NaN],'YData',[P(:,2); NaN],'ZData',[P(:,3); NaN], ...
                        'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','EdgeAlpha',trial_alpha,'LineWidth',0.5)
                else
                    patch('XData',[P(:,1); NaN],'YData',[P(:,2); NaN], ...
                        'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','EdgeAlpha',trial_alpha,'LineWidth',0.5)
                end
            end
            %invisible lines, just so the legend has solid colours
            hgrp = gobjects(1,2);
            for g = 1:2
                hgrp(g) = plot(nan,nan,'Color',grpcol(g,:),'LineWidth',2);
            end
            hold off
            legend(hgrp,{['no stim, ' num2str(sum(with_stim{prb} == false)) ' events'],['with stim, ' num2str(sum(with_stim{prb})) ' events']},'Location','best')
            title([char(signames(j)) ', ' num2str(nch) ' channels, every event'],'Interpreter','none')

            %bottom panel: the means
            nexttile(nband+j)
            hold on
            hmean = gobjects(0);
            labs = {};
            for g = 1:2 %event-start-aligned mean of each group
                evs = find(with_stim{prb} == (g == 2))';
                if isempty(evs) == true
                    continue
                end
                M = zeros(nbp,nch); %[bins x channels]
                for ev = evs
                    i1 = find(binstart{prb}{ev} >= ev_times(ev,1),1); %first bin of the event
                    M = M + double(envz{prb}{j}{ev}(i1:i1+nbp-1,:)); %snippets reach 1 s past the event, so this always fits
                end
                P = (M/numel(evs) - mu)*coeff; %[bins x pca_dims]
                gc = grpcol(g,:);
                fadecol = gc + fade_start*(1 - gc);
                tracecol = fadecol + linspace(0,1,nbp)'*(gc - fadecol);
                if pca_dims == 3
                    patch('XData',[P(:,1); NaN],'YData',[P(:,2); NaN],'ZData',[P(:,3); NaN],'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','LineWidth',2)
                    plot3(P(1,1),P(1,2),P(1,3),'o','Color',gc,'MarkerFaceColor',gc)
                    plot3(P(end,1),P(end,2),P(end,3),'s','Color',gc,'MarkerFaceColor',gc,'MarkerSize',8)
                else
                    patch('XData',[P(:,1); NaN],'YData',[P(:,2); NaN],'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','LineWidth',2)
                    plot(P(1,1),P(1,2),'o','Color',gc,'MarkerFaceColor',gc)
                    plot(P(end,1),P(end,2),'s','Color',gc,'MarkerFaceColor',gc,'MarkerSize',8)
                end
                hmean(end+1) = plot(nan,nan,'Color',gc,'LineWidth',2); %#ok<SAGROW>
                labs{end+1} = [group_names{g} ' events, first ' num2str(ev_plot_ms) ' ms (n=' num2str(numel(evs)) ')']; %#ok<SAGROW>
            end
            for cnd = 1:nstimcond %stim-aligned mean of each stim condition
                idx = find(strcmp(stim_condition,stim_conds{cnd}) & isnan(stim_ev{prb}) == false)';
                if isempty(idx) == true
                    continue
                end
                M = zeros(nbs,nch);
                for s = idx
                    ev = stim_ev{prb}(s);
                    i1 = find(binstart{prb}{ev} >= stim_onset(s) + pre_win(1),1); %first bin of the pre window
                    M = M + double(envz{prb}{j}{ev}(i1:i1+nbs-1,:)); %stim_ev already checked the snippet covers this
                end
                P = (M/numel(idx) - mu)*coeff;
                gc = condcol(cnd,:);
                fadecol = gc + fade_start*(1 - gc);
                tracecol = fadecol + linspace(0,1,nbs)'*(gc - fadecol);
                if pca_dims == 3
                    patch('XData',[P(:,1); NaN],'YData',[P(:,2); NaN],'ZData',[P(:,3); NaN],'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','LineWidth',2)
                    plot3(P(onset_bs,1),P(onset_bs,2),P(onset_bs,3),'o','Color',gc,'MarkerFaceColor',gc)
                    plot3(P(end,1),P(end,2),P(end,3),'s','Color',gc,'MarkerFaceColor',gc,'MarkerSize',8)
                else
                    patch('XData',[P(:,1); NaN],'YData',[P(:,2); NaN],'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','LineWidth',2)
                    plot(P(onset_bs,1),P(onset_bs,2),'o','Color',gc,'MarkerFaceColor',gc)
                    plot(P(end,1),P(end,2),'s','Color',gc,'MarkerFaceColor',gc,'MarkerSize',8)
                end
                hmean(end+1) = plot(nan,nan,'Color',gc,'LineWidth',2); %#ok<SAGROW>
                labs{end+1} = [stim_conds{cnd} ', ' num2str(pre_win(1)) ' to ' num2str(post_win(2)) ' ms re: stim (n=' num2str(numel(idx)) ')']; %#ok<SAGROW>
            end
            hold off
            if isempty(hmean) == false
                legend(hmean,labs,'Location','best','Interpreter','none','FontSize',7)
            end
            title([char(signames(j)) ', means'],'Interpreter','none')

            %same axis labels on both panels
            for panel = [j nband+j]
                nexttile(panel)
                xlabel(['PC1 (' num2str(pcexpl(1),'%.1f') '%)'])
                ylabel(['PC2 (' num2str(pcexpl(2),'%.1f') '%)'])
                if pca_dims == 3
                    zlabel(['PC3 (' num2str(pcexpl(3),'%.1f') '%)'])
                    view(3)
                    grid on
                end
            end
        end
        sgtitle({[animal ' - probe ' num2str(probes(prb)) ': events in PC space, grey = no stim, red = with stim'], ...
            'direction: faded = early, full colour = late, dot = event start or stim onset, square = end'},'Interpreter','none')
    end
end


%% Data Visualization: Stim-aligned trajectories in PC space, per recording
%one figure per probe and recording with stims, one column per band. Every stim whose onset lands inside an event is a line
%from align_win(1) to align_win(2) around the onset: grey = before the stim, red = after. Thick line = their mean, dot = stim onset, square = end.
%align_shift = 1 moves every line so its onset point sits at 0, align_pcfit picks one PCA for all recordings or one per recording.
%every recording's panel of a band gets the same axis limits (the widest any recording needs), so the space they take up compares directly
%the PCA is rebuilt here, since the PCA plot section above can be off
alignTraj = cell(length(probes),nfiles); %{probe x recording}{band} [stims x bins x numel(align_pcs)] each stim's line through the PCs in align_pcs
alignStim = cell(length(probes),nfiles); %{probe x recording} [stims x 1] which stims (index into stim_onset) those lines are
if plot_stimalign == 1
    disp('PCA: Lining up stims in PC space')
    nba = round((align_win(2)-align_win(1))/bin_ms); %bins in each stim's line
    onset_ba = round(-align_win(1)/bin_ms) + 1; %bin where the stim comes on (within one bin_ms)
    shifttxt = {'not shifted','shifted so every onset sits at 0'};
    fittxt = {'one PCA for all recordings','one PCA per recording'};
    npc = numel(align_pcs); %how many PCs get plotted, 2 or 3
    for prb = 1:length(probes)
        ev_times = oscillations.prbEvent_timestamps{prb}; %[events x 2], start and stop samples
        ev_ext = oscillations.prbEvent_ext_timestamps{prb}; %[events x 2], what the snippets cover
        ev_file = oscillations.prbEvent_file_index{prb};
        nev = size(ev_times,1);

        %which event does each stim land in? stays NaN if none, or if that event's snippet doesn't cover align_win + edge_ms
        aev = nan(nstim,1);
        for s = 1:nstim
            ev = find(ev_file == stim_file(s) & ev_times(:,1) <= stim_onset(s) & ev_times(:,2) >= stim_onset(s) ...
                & ev_ext(:,1) <= stim_onset(s) + align_win(1) - edge_ms & ev_ext(:,2) >= stim_onset(s) + align_win(2) - 1 + edge_ms, 1);
            if isempty(ev) == false
                aev(s) = ev;
            end
        end

        %one PCA per band on the bins inside every event, all recordings together
        if align_pcfit == 0
            apc = cell(1,nband); %{band} struct: coeff [channels x npc] (only the align_pcs columns), mu [1 x channels], expl [PCs x 1]
            for j = 1:nband
                Xall = [];
                for ev = 1:nev
                    Xall = [Xall; double(envz{prb}{j}{ev}(evmask{prb}{ev},:))]; %#ok<AGROW>
                end
                [coeff,~,latent,~,~,mu] = pca(Xall,'NumComponents',max(align_pcs));
                apc{j} = struct('coeff',coeff(:,align_pcs),'mu',mu,'expl',100*latent/sum(latent)); %keep just the PCs you asked for
            end
            clear Xall
        end

        ax_align = gobjects(nfiles,nband); %[recordings x bands] each panel, so they can all get the same axis limits at the end
        for file = 1:nfiles
            idx = find(stim_file == file & isnan(aev) == false)'; %this recording's stims that land inside an event
            if isempty(idx) == true
                if any(stim_file == file) == true %spon recordings have no stims at all, so they don't get a note
                    run_notes(end+1) = "Probe " + probes(prb) + ", recording " + file + ": no stim-aligned PCA plot, no stim landed inside an event"; %#ok<SAGROW>
                end
                continue
            end

            %or one PCA per band on just this recording's events
            if align_pcfit == 1
                apc = cell(1,nband);
                for j = 1:nband
                    Xall = [];
                    for ev = find(ev_file == file)'
                        Xall = [Xall; double(envz{prb}{j}{ev}(evmask{prb}{ev},:))]; %#ok<AGROW>
                    end
                    [coeff,~,latent,~,~,mu] = pca(Xall,'NumComponents',max(align_pcs));
                    apc{j} = struct('coeff',coeff(:,align_pcs),'mu',mu,'expl',100*latent/sum(latent)); %keep just the PCs you asked for
                end
                clear Xall
            end

            alignStim{prb,file} = idx(:);
            alignTraj{prb,file} = cell(1,nband);
            figure('Name',[animal '-probe' num2str(probes(prb)) '-rec' num2str(file) '-stimalignPCA'],'Color','w');
            set(gcf, 'Position', get(0, 'Screensize'));
            tiledlayout(1,nband) %one column per band

            for j = 1:nband
                %every stim's line through PC space
                P = zeros(numel(idx),nba,npc); %[stims x bins x npc]
                for k = 1:numel(idx)
                    s = idx(k);
                    ev = aev(s);
                    i1 = find(binstart{prb}{ev} >= stim_onset(s) + align_win(1),1); %first bin of the window, the coverage check above makes sure the rest fits
                    Pk = (double(envz{prb}{j}{ev}(i1:i1+nba-1,:)) - apc{j}.mu)*apc{j}.coeff; %[bins x npc]
                    if align_shift == 1
                        Pk = Pk - Pk(onset_ba,:); %onset point to 0
                    end
                    P(k,:,:) = Pk;
                end
                alignTraj{prb,file}{j} = P;
                Pm = reshape(mean(P,1),nba,npc); %[bins x npc] mean line

                %thin see-through lines for each stim, grey up to the onset and red from it (both include the onset bin so there's no gap)
                ax_align(file,j) = nexttile(j);
                hold on
                for k = randperm(numel(idx)) %random drawing order so no stim always sits on top
                    Pk = reshape(P(k,:,:),nba,npc);
                    if npc == 3
                        plot3(Pk(1:onset_ba,1),Pk(1:onset_ba,2),Pk(1:onset_ba,3),'Color',[grpcol(1,:) align_alpha],'LineWidth',align_width)
                        plot3(Pk(onset_ba:end,1),Pk(onset_ba:end,2),Pk(onset_ba:end,3),'Color',[grpcol(2,:) align_alpha],'LineWidth',align_width)
                    else
                        plot(Pk(1:onset_ba,1),Pk(1:onset_ba,2),'Color',[grpcol(1,:) align_alpha],'LineWidth',align_width)
                        plot(Pk(onset_ba:end,1),Pk(onset_ba:end,2),'Color',[grpcol(2,:) align_alpha],'LineWidth',align_width)
                    end
                end
                %thick mean on top, dot = onset, square = end
                if npc == 3
                    hpre = plot3(Pm(1:onset_ba,1),Pm(1:onset_ba,2),Pm(1:onset_ba,3),'Color',grpcol(1,:),'LineWidth',2.5);
                    hpost = plot3(Pm(onset_ba:end,1),Pm(onset_ba:end,2),Pm(onset_ba:end,3),'Color',grpcol(2,:),'LineWidth',2.5);
                    plot3(Pm(onset_ba,1),Pm(onset_ba,2),Pm(onset_ba,3),'o','Color','k','MarkerFaceColor','k')
                    plot3(Pm(end,1),Pm(end,2),Pm(end,3),'s','Color',grpcol(2,:),'MarkerFaceColor',grpcol(2,:),'MarkerSize',8)
                else
                    hpre = plot(Pm(1:onset_ba,1),Pm(1:onset_ba,2),'Color',grpcol(1,:),'LineWidth',2.5);
                    hpost = plot(Pm(onset_ba:end,1),Pm(onset_ba:end,2),'Color',grpcol(2,:),'LineWidth',2.5);
                    plot(Pm(onset_ba,1),Pm(onset_ba,2),'o','Color','k','MarkerFaceColor','k')
                    plot(Pm(end,1),Pm(end,2),'s','Color',grpcol(2,:),'MarkerFaceColor',grpcol(2,:),'MarkerSize',8)
                end
                hold off
                legend([hpre hpost],{[num2str(align_win(1)) ' to 0 ms (before stim)'],['0 to ' num2str(align_win(2)) ' ms (after stim)']},'Location','best')
                title([char(signames(j)) ', ' num2str(numel(idx)) ' stims'],'Interpreter','none')
                xlabel(['PC' num2str(align_pcs(1)) ' (' num2str(apc{j}.expl(align_pcs(1)),'%.1f') '%)'])
                ylabel(['PC' num2str(align_pcs(2)) ' (' num2str(apc{j}.expl(align_pcs(2)),'%.1f') '%)'])
                if npc == 3
                    zlabel(['PC' num2str(align_pcs(3)) ' (' num2str(apc{j}.expl(align_pcs(3)),'%.1f') '%)'])
                    view(3)
                    grid on
                end
            end
            sgtitle({[animal ' - probe ' num2str(probes(prb)) ' - recording ' num2str(file) ' (' file_condition{file} '): stims inside events, lined up on stim onset'], ...
                ['grey = before the stim, red = after, thick = mean, dot = onset, square = end. ' shifttxt{align_shift+1} ', ' fittxt{align_pcfit+1}]},'Interpreter','none')
        end

        %same axis limits on every recording's panel of a band, so you can compare how much space they take up.
        %limits = smallest and biggest value any recording needs on that PC, plus 5% so nothing sits on the edge
        for j = 1:nband
            lo = inf(1,npc);
            hi = -inf(1,npc);
            for file = 1:nfiles
                if isempty(alignTraj{prb,file}) == false
                    lo = min(lo, reshape(min(alignTraj{prb,file}{j},[],[1 2]),1,npc)); %smallest value on each PC
                    hi = max(hi, reshape(max(alignTraj{prb,file}{j},[],[1 2]),1,npc)); %biggest value on each PC
                end
            end
            pad = 0.05*(hi - lo);
            for file = find(isgraphics(ax_align(:,j)))'
                xlim(ax_align(file,j),[lo(1)-pad(1) hi(1)+pad(1)])
                ylim(ax_align(file,j),[lo(2)-pad(2) hi(2)+pad(2)])
                if npc == 3
                    zlim(ax_align(file,j),[lo(3)-pad(3) hi(3)+pad(3)])
                end
            end
        end
    end
end


%% Data Visualization: Scroll through events (left/right = event, up/down = colors)
%only hide events from the viewer's copies, the real oscillations and prbsnippets stay whole
osc_view = oscillations;
snip_view = prbsnippets;
if view_which > 0
    for prb = 1:length(probes)
        keep = with_stim{prb} == (view_which == 1);
        osc_view.prbEvent_file_index{prb} = osc_view.prbEvent_file_index{prb}(keep);
        osc_view.prbEvent_ext_timestamps{prb} = osc_view.prbEvent_ext_timestamps{prb}(keep,:);
        snip_view.event_data{prb} = snip_view.event_data{prb}(keep);
    end
end
oscillation_plotEvents(ProbeInfo,osc_view,snip_view)


%% Data Visualization: Scroll through bouts with their PR over time and manifold
%like the viewer above, plus each bout's sliding PR (top) and a 3-PC manifold of it (right)
%every probe's bouts in one time-ordered list. v shows probe 1, probe 2 or all of them, over exactly the bout's time
%keys: left/right = bout, v = which probes, up/down = colours (shank/depth), b = band, p = own/shared PCs, e = LFP/envelope traces
%the key and drawing code are local functions at the very end of this script
if view_pr == 1
    %which bouts? same rule as the viewer above
    vflat = []; %[bouts x 2] probe index, event index
    for prb = 1:length(probes)
        evs = (1:numel(with_stim{prb}))';
        if view_which > 0
            evs = evs(with_stim{prb} == (view_which == 1));
        end
        vflat = [vflat; [repmat(prb,numel(evs),1) evs]]; %#ok<AGROW>
    end
    %put every probe's bouts in time order: recording first, then bout start
    vkey = zeros(size(vflat,1),2); %[bouts x 2] recording, start sample
    for b = 1:size(vflat,1)
        vkey(b,:) = [oscillations.prbEvent_file_index{vflat(b,1)}(vflat(b,2)) oscillations.prbEvent_timestamps{vflat(b,1)}(vflat(b,2),1)];
    end
    [~,vord] = sortrows(vkey);
    vflat = vflat(vord,:);

    if isempty(vflat) == true
        run_notes(end+1) = "No bout viewer: no events to show";
    else
        %shared PCs: one PCA per probe and band on the bins inside every event (rebuilt here, the PCA plot section can be off)
        vpc = cell(length(probes),nband); %{probe x band} struct: coeff [channels x 3], mu [1 x channels], expl [PCs x 1]
        for prb = 1:length(probes)
            for j = 1:nband
                Xall = [];
                for ev = 1:numel(envz{prb}{j})
                    Xall = [Xall; double(envz{prb}{j}{ev}(evmask{prb}{ev},:))]; %#ok<AGROW>
                end
                [coeff,~,latent,~,~,mu] = pca(Xall,'NumComponents',3);
                vpc{prb,j} = struct('coeff',coeff,'mu',mu,'expl',100*latent/sum(latent));
            end
        end
        clear Xall

        %which shank and how deep is each channel? same lookup as oscillation_plotEvents.m (ProbeMaps column = shank, row after Ch_Remove = depth)
        vchans = cell(1,length(probes)); %{probe} rows of the snippets that belong to this probe
        vshank = cell(1,length(probes)); %{probe} [channels x 1] shank of each channel
        vdepth = cell(1,length(probes)); %{probe} [channels x 1] position along its shank
        for prb = 1:length(probes)
            vchans{prb} = find(oscillations.chan_probe == prb);
            ids = oscillations.chanals(vchans{prb});
            pm = ProbeInfo.ProbeMaps{probes(prb)};
            chrem = [];
            if numel(ProbeInfo.Ch_Remove) >= probes(prb)
                chrem = ProbeInfo.Ch_Remove{probes(prb)};
            end
            vshank{prb} = zeros(numel(ids),1);
            vdepth{prb} = zeros(numel(ids),1);
            for i = 1:numel(ids)
                [~,col] = find(pm == ids(i));
                vshank{prb}(i) = col(1);
                kept = pm(:,col(1));
                kept = kept(ismember(kept,chrem) == false & kept ~= 0); %this shank's channels, physical order, removed ones left out
                vdepth{prb}(i) = find(kept == ids(i),1);
            end
        end

        %everything the key and drawing functions need goes into the figure
        ud = struct();
        ud.flat = vflat;
        ud.cur = 1; %which bout
        ud.band = 1; %which band for the manifold and envelope traces
        ud.pcmode = 'own'; %'own' = this bout's PCs, 'shared' = PCs of all bouts
        ud.tracemode = 'lfp'; %'lfp' = 2-80 Hz LFP, 'env' = envelope of the selected band
        ud.cmapIdx = 1; %1 = shank colours, 2 = depth colours
        ud.nprb = length(probes);
        ud.show = 0; %which probes to draw: 0 = all of them, k = only probe k (v key cycles it)
        if ud.nprb == 1
            ud.show = 1;
        end
        ud.hbert_sos = hbert_sos; %band filters, for the other probes' envelopes (envz only has each bout's own probe)
        ud.hbert_g = hbert_g;
        ud.sdch = sdch; %SDs for the z-score, so those envelopes are in the same z-units as envz
        ud.envz = envz;
        ud.binstart = binstart;
        ud.evmask = evmask;
        ud.snips = prbsnippets.event_data;
        ud.chans = vchans;
        ud.shank = vshank;
        ud.depth = vdepth;
        ud.pc = vpc;
        ud.ev_times = oscillations.prbEvent_timestamps;
        ud.ev_ext = oscillations.prbEvent_ext_timestamps;
        ud.ev_file = oscillations.prbEvent_file_index;
        ud.stim_times = oscillations.stim_times;
        ud.with_stim = with_stim;
        ud.eventTable = eventTable;
        ud.probes = probes;
        ud.signames = signames;
        ud.bandtags = bandtags;
        ud.nband = nband;
        ud.bin_ms = bin_ms;
        ud.slide_ms = slide_ms;
        ud.slide_step_ms = slide_step_ms;
        ud.edge_ms = edge_ms;
        ud.stim_reach = stim_reach;
        ud.pr_calc = pr_calc;
        ud.fade_start = fade_start;

        fig = figure('Name',[animal ' bout viewer: PR and manifold'],'Color','w','KeyPressFcn',@oscpr_keypress);
        set(fig, 'Position', get(0, 'Screensize'));
        ud.axPR = axes(fig,'Position',[0.05 0.70 0.58 0.18]); %PR over time, every shown probe
        %one probe shown: traces and manifold fill the space under the PR
        ud.axTr = axes(fig,'Position',[0.05 0.06 0.58 0.60]); %stacked traces
        ud.axMf = axes(fig,'Position',[0.69 0.15 0.29 0.65]); %manifold
        view(ud.axMf,3) %set once, so a rotation you make stays when you change bout
        %all probes shown: one row per probe under the PR, probe 1 on top
        rowh = 0.62/ud.nprb; %height of each row
        ud.axTrN = gobjects(1,ud.nprb);
        ud.axMfN = gobjects(1,ud.nprb);
        for q = 1:ud.nprb
            y0 = 0.04 + (ud.nprb-q)*rowh;
            ud.axTrN(q) = axes(fig,'Position',[0.05 y0 0.58 rowh-0.04]);
            ud.axMfN(q) = axes(fig,'Position',[0.69 y0 0.29 rowh-0.04]);
            view(ud.axMfN(q),3)
        end
        setappdata(fig,'ud',ud);
        oscpr_redraw(fig);
    end
end


%% Post-run report
disp('---------- Oscillation PR report ----------')
disp(fileReport)
for prb = 1:length(probes)
    disp(['Probe ' num2str(probes(prb)) ', ' num2str(nch_prb(prb)) ' channels'])
    disp(stimEventReport{prb})
    disp(stimSortSummary{prb})
    disp(prepostSummary{prb})
end
disp('PR: per-window is capped by how few independent points a window holds (worst in narrow, slow bands). Pooled isn''t, but loud windows count for more. Equal-weight pooled counts every window the same. n_indep should be well above the channel count.')
for i = 1:numel(run_notes)
    disp(run_notes(i))
end


%Now the local functions for the bout viewer. MATLAB wants these at the very end of the script

function oscpr_keypress(fig,evt)
%keys for the bout viewer: left/right = bout, up/down = colours, b = band, p = own/shared PCs, e = LFP/envelope traces, v = which probes
ud = getappdata(fig,'ud');
switch evt.Key
    case 'v'
        if ud.nprb > 1
            ud.show = mod(ud.show + 1, ud.nprb + 1); %probe 1, probe 2, ..., then 0 = all of them
        end
    case 'rightarrow'
        ud.cur = min(ud.cur + 1, size(ud.flat,1));
    case 'leftarrow'
        ud.cur = max(ud.cur - 1, 1);
    case {'uparrow','downarrow'}
        ud.cmapIdx = 3 - ud.cmapIdx; %only two schemes, shank and depth
    case 'b'
        ud.band = mod(ud.band, ud.nband) + 1;
    case 'p'
        if strcmp(ud.pcmode,'own') == true
            ud.pcmode = 'shared';
        else
            ud.pcmode = 'own';
        end
    case 'e'
        if strcmp(ud.tracemode,'lfp') == true
            ud.tracemode = 'env';
        else
            ud.tracemode = 'lfp';
        end
    otherwise
        return
end
setappdata(fig,'ud',ud);
oscpr_redraw(fig);
end


function oscpr_redraw(fig)
%draws the current bout: sliding PR of every shown probe on top, then each shown probe's stacked traces and manifold
ud = getappdata(fig,'ud');
prb = ud.flat(ud.cur,1); %the probe that found this bout
ev = ud.flat(ud.cur,2);
j = ud.band;
ext = ud.ev_ext{prb}(ev,:); %snippet start and stop, recording samples (1 kHz, so ms)
evt = ud.ev_times{prb}(ev,:); %bout start and stop
file = ud.ev_file{prb}(ev);
sb = ud.binstart{prb}{ev}; %[bins x 1] recording sample where each bin starts
bc = sb + (ud.bin_ms-1)/2; %bin centres
okbin = sb >= ext(1) + ud.edge_ms & sb + ud.bin_ms - 1 <= ext(2) - ud.edge_ms; %bins clear of the snippet edges
inev = ud.evmask{prb}{ev}; %bins inside the bout
bandcol = lines(ud.nband);
if ud.with_stim{prb}(ev) == true
    gc = [0.85 0.15 0.15]; %with stim = red
    stimtxt = 'with stim';
else
    gc = [0 0 0]; %no stim = black
    stimtxt = 'no stim';
end

%stims: red lines for the ones inside the snippet, grey shading up to stim_reach after any stim (same as oscillation_plotEvents.m)
st = ud.stim_times{file};
st = st(:);
st_in = st(st >= ext(1) & st <= ext(2));
abs_t = ext(1):ext(2);
near = false(size(abs_t));
for s = st(st >= ext(1) - ud.stim_reach & st <= ext(2))'
    near((max(ext(1),s):min(ext(2),s + ud.stim_reach)) - ext(1) + 1) = true;
end
d = diff([0 near 0]);
run_starts = abs_t(d(1:end-1) == 1); %first sample of each shaded run
run_ends = abs_t(find(d == -1) - 1);

%which probes to draw
if ud.show == 0
    shown = 1:ud.nprb;
else
    shown = ud.show;
end

%every shown probe's binned z-scored envelopes over EXACTLY this bout's snippet, {probe x band} [bins x channels]
%the bout's own probe already has them in envz. The snippet holds every analyzed channel, so the other probes get
%the same recipe as the envelope section, from the same samples, and the same SDs as their own envz
E = cell(ud.nprb,ud.nband);
nbin = numel(sb);
for q = shown
    if q == prb
        for jj = 1:ud.nband
            E{q,jj} = double(ud.envz{prb}{jj}{ev});
        end
    else
        snip = ud.snips{prb}{ev}(ud.chans{q},:)'; %[snippet samples x channels] probe q, same samples as the bout
        for jj = 1:ud.nband
            env = abs(hilbert(filtfilt(ud.hbert_sos{jj}, ud.hbert_g(jj), snip))); %[snippet samples x channels] amplitude envelope, uV
            env = env(1:nbin*ud.bin_ms,:); %drop the leftover samples that don't fill a bin
            E{q,jj} = reshape(mean(reshape(env,ud.bin_ms,nbin,[]),1),nbin,[]) ./ ud.sdch{q}(jj,:); %[bins x channels] z-units, bins line up with the bout's
        end
    end
end

%---- top: sliding per-window PR, each probe on its own channels. Bout's probe solid, other probes dashed ----
ax = ud.axPR;
cla(ax);
hold(ax,'on');
starts = (ext(1) + ud.edge_ms):ud.slide_step_ms:(ext(2) - ud.edge_ms - ud.slide_ms + 1); %window starts
prts = nan(ud.nprb,ud.nband,numel(starts)); %[probes x bands x windows]
for q = shown
    for jj = 1:ud.nband
        for t = 1:numel(starts)
            inw = sb >= starts(t) & sb + ud.bin_ms - 1 <= starts(t) + ud.slide_ms - 1; %bins in this window
            prts(q,jj,t) = ud.pr_calc(cov(E{q,jj}(inw,:)));
        end
    end
end
yhi = max([prts(:); 2])*1.1; %max skips the NaNs of probes not shown
patch(ax,[evt(1) evt(2) evt(2) evt(1)],[1 1 yhi yhi],[1 0.93 0.7],'EdgeColor','none','HandleVisibility','off'); %the bout
for r = 1:numel(run_starts)
    patch(ax,[run_starts(r) run_ends(r) run_ends(r) run_starts(r)],[1 1 yhi yhi],[0.7 0.7 0.7],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off');
end
for q = shown
    if q == prb
        ls = '-';
        probetxt = " (this bout)";
    else
        ls = '--';
        probetxt = " (dashed)";
    end
    for jj = 1:ud.nband
        if jj == j
            lw = 2.5; %selected band thick
        else
            lw = 1;
        end
        plot(ax,starts + ud.slide_ms/2,reshape(prts(q,jj,:),1,[]),ls,'Color',bandcol(jj,:),'LineWidth',lw, ...
            'DisplayName',ud.bandtags(jj) + " Hz, probe " + ud.probes(q) + probetxt);
    end
end
for s = st_in'
    xline(ax,s,'Color','r','LineWidth',1.5,'HandleVisibility','off');
end
hold(ax,'off');
xlim(ax,[ext(1) ext(2)]);
ylim(ax,[1 yhi]);
ylabel(ax,'PR');
legend(ax,'Location','northeastoutside','FontSize',7);
set(ax,'Position',[0.05 0.70 0.58 0.18]); %the outside legend shrinks the axes, so put it back to the traces' width and the times line up
title(ax,sprintf('per-window PR per probe, %d ms window, %d ms steps (thick = selected band, solid = this bout''s probe, dashed = other probe, yellow = bout, grey = after a stim)',ud.slide_ms,ud.slide_step_ms),'FontWeight','normal');

%---- traces and manifold of each shown probe ----
%one probe shown = the big pair of axes, more = one row per probe. The pair not in use is wiped and hidden
if numel(shown) == 1
    axTr = ud.axTr;
    axMf = ud.axMf;
    unused = [ud.axTrN ud.axMfN];
else
    axTr = ud.axTrN(shown);
    axMf = ud.axMfN(shown);
    unused = [ud.axTr ud.axMf];
end
for a = unused
    if strcmp(a.Visible,'on') == true
        cla(a,'reset'); %reset also clears titles and labels, which plain cla leaves showing
        a.Visible = 'off';
    end
end
%everything about the window that every probe's panels share
w = struct('prb',prb,'ev',ev,'j',j,'ext',ext,'evt',evt,'file',file,'sb',sb,'bc',bc,'okbin',okbin,'inev',inev, ...
    'abs_t',abs_t,'st_in',st_in,'run_starts',run_starts,'run_ends',run_ends,'gc',gc,'compact',numel(shown) > 1);
for k = 1:numel(shown)
    oscpr_drawprobe(ud, shown(k), E(shown(k),:), axTr(k), axMf(k), w);
end

%---- what's on screen ----
if ud.cmapIdx == 1
    schemetxt = 'shank';
else
    schemetxt = 'depth';
end
if numel(shown) == 1
    showtxt = ['probe ' num2str(ud.probes(shown))];
else
    showtxt = 'all probes';
end
et = ud.eventTable{prb};
prtxt = strings(1,ud.nband);
for jj = 1:ud.nband
    prtxt(jj) = ud.bandtags(jj) + " Hz: " + num2str(et.(char("pr_" + ud.bandtags(jj)))(ev),'%.2f');
end
sgtitle(fig,{sprintf('Bout found on probe %d - Recording %d (%s): bout %d/%d on that probe (overall %d/%d), %s, %d ms', ...
    ud.probes(prb), file, string(et.condition(ev)), ev, height(et), ud.cur, size(ud.flat,1), stimtxt, round(et.dur_ms(ev))), ...
    "bout PR on probe " + ud.probes(prb) + "  " + strjoin(prtxt,",  "), ...
    ['showing: ' showtxt ' | band: ' char(ud.bandtags(j)) ' Hz | PCs: ' ud.pcmode ' | traces: ' ud.tracemode ' | colours: ' schemetxt ...
    '      keys: left/right bout, v probes, up/down colours, b band, p PCs, e traces']},'Interpreter','none','FontSize',10);
end


function oscpr_drawprobe(ud, q, Eq, axTr, axMf, w)
%draws one probe's stacked traces and manifold over the current bout's window
%q = which probe to draw, Eq = its {1 x band} [bins x channels] z-scored envelopes over that window, w = the window info from oscpr_redraw
%every probe uses the bout's own snippet, so the samples are exactly the same on every probe
if q == w.prb
    rowtxt = ['probe ' num2str(ud.probes(q)) ' (this bout)'];
else
    rowtxt = ['probe ' num2str(ud.probes(q)) ' (same time)'];
end
%axes that were hidden got reset, so switch them back on (and the manifold back to 3D)
if strcmp(axTr.Visible,'off') == true
    axTr.Visible = 'on';
end
if strcmp(axMf.Visible,'off') == true
    axMf.Visible = 'on';
    view(axMf,3)
end

%---- left: stacked traces ----
ax = axTr;
cla(ax);
hold(ax,'on');
nch = numel(ud.chans{q});
if ud.cmapIdx == 1 %shank colours, natural channel order (already grouped by shank)
    cols = lines(max(ud.shank{q}));
    chcol = cols(ud.shank{q},:);
    order = (1:nch)';
    schemetxt = 'shank';
else %depth colours, grouped by depth first and shank second
    cols = copper(max(ud.depth{q}));
    chcol = cols(ud.depth{q},:);
    [~,order] = sortrows([ud.depth{q} ud.shank{q} (1:nch)']);
    schemetxt = 'depth';
end
if strcmp(ud.tracemode,'lfp') == true
    Y = ud.snips{w.prb}{w.ev}(ud.chans{q},:); %[channels x samples] 2-80 Hz LFP of probe q, from the bout's own snippet, uV
    tt = w.abs_t;
    step = 250; %uV between channels, same as oscillation_plotEvents.m
    tracetxt = '2-80 Hz LFP (uV)';
else
    Y = Eq{w.j}'; %[channels x bins] z-scored envelope
    tt = w.bc';
    step = 4; %z-units between channels
    tracetxt = [char(ud.bandtags(w.j)) ' Hz envelope (z)'];
end
ylo = min(Y(:)) - step*(nch-1);
ytop = max(Y(:));
patch(ax,[w.evt(1) w.evt(2) w.evt(2) w.evt(1)],[ylo ylo ytop ytop],[1 0.93 0.7],'EdgeColor','none'); %the bout, yellow on every probe
%light blue = this probe's own bouts in the window, so you can see if it found one here too
if q ~= w.prb
    for o = find(ud.ev_file{q} == w.file & ud.ev_times{q}(:,1) <= w.ext(2) & ud.ev_times{q}(:,2) >= w.ext(1))'
        o1 = max(w.ext(1),ud.ev_times{q}(o,1));
        o2 = min(w.ext(2),ud.ev_times{q}(o,2));
        patch(ax,[o1 o2 o2 o1],[ylo ylo ytop ytop],[0.6 0.8 1],'FaceAlpha',0.4,'EdgeColor','none');
    end
end
for r = 1:numel(w.run_starts)
    patch(ax,[w.run_starts(r) w.run_ends(r) w.run_ends(r) w.run_starts(r)],[ylo ylo ytop ytop],[0.7 0.7 0.7],'FaceAlpha',0.3,'EdgeColor','none');
end
for slot = 1:nch
    i = order(slot);
    plot(ax,tt,Y(i,:) - step*(slot-1),'Color',chcol(i,:));
end
for s = w.st_in'
    xline(ax,s,'Color','r','LineWidth',1.5);
end
hold(ax,'off');
xlim(ax,[w.ext(1) w.ext(2)]);
ylim(ax,[ylo ytop]);
set(ax,'YTick',[]);
xlabel(ax,'Time in Recording (ms)');
title(ax,[rowtxt ', traces: ' tracetxt ', colours: ' schemetxt ' (yellow = bout, blue = this probe''s own bouts)'],'FontWeight','normal');

%---- right: the bout's stretch as a line through the top 3 PCs ----
ax = axMf;
cla(ax);
hold(ax,'on');
X = Eq{w.j}; %[bins x channels]
if strcmp(ud.pcmode,'own') == true
    [coeff,~,latent,~,~,mu] = pca(X(w.inev,:),'NumComponents',3); %PCs of probe q during the bout
    expl = 100*latent/sum(latent);
    pctxt = 'this bout''s PCs';
else
    coeff = ud.pc{q,w.j}.coeff; %PCs of all of probe q's bouts
    mu = ud.pc{q,w.j}.mu;
    expl = ud.pc{q,w.j}.expl;
    pctxt = 'shared PCs (all bouts)';
end
P = (X - mu)*coeff; %[bins x 3]
ib = find(w.inev); %bout bins (the same bins on every probe)
pre = find(w.okbin & w.sb < w.evt(1)); %before the bout, edges left out
post = find(w.okbin & w.inev == false & w.sb > w.evt(1)); %after the bout
if isempty(pre) == false
    pre = [pre; ib(1)]; %join up with the bout
    plot3(ax,P(pre,1),P(pre,2),P(pre,3),'-','Color',[0.75 0.75 0.75],'LineWidth',1);
end
if isempty(post) == false
    post = [ib(end); post];
    plot3(ax,P(post,1),P(post,2),P(post,3),'-','Color',[0.75 0.75 0.75],'LineWidth',1);
end
fadecol = w.gc + ud.fade_start*(1 - w.gc); %washed out version of the bout colour, its start
tracecol = fadecol + linspace(0,1,numel(ib))'*(w.gc - fadecol); %[bout bins x RGB] faded to full colour
%a patch with no face and interpolated edge colours is a line whose colour changes along it. The NaN stops it closing into a loop
patch(ax,'XData',[P(ib,1); NaN],'YData',[P(ib,2); NaN],'ZData',[P(ib,3); NaN], ...
    'FaceVertexCData',[tracecol; 1 1 1],'FaceColor','none','EdgeColor','interp','LineWidth',2.5);
plot3(ax,P(ib(1),1),P(ib(1),2),P(ib(1),3),'o','Color',w.gc,'MarkerFaceColor',w.gc,'MarkerSize',8);
plot3(ax,P(ib(end),1),P(ib(end),2),P(ib(end),3),'s','Color',w.gc,'MarkerFaceColor',w.gc,'MarkerSize',9);
for s = w.st_in'
    [~,b] = min(abs(w.bc - s)); %bin nearest the stim
    if w.okbin(b) == true
        plot3(ax,P(b,1),P(b,2),P(b,3),'p','Color','r','MarkerFaceColor','r','MarkerSize',14);
    end
end
hold(ax,'off');
grid(ax,'on');
xlabel(ax,['PC1 (' num2str(expl(1),'%.1f') '%)']);
ylabel(ax,['PC2 (' num2str(expl(2),'%.1f') '%)']);
zlabel(ax,['PC3 (' num2str(expl(3),'%.1f') '%)']);
if w.compact == true %rows are short when every probe is shown, so keep the title to one line
    title(ax,[rowtxt ', ' char(ud.signames(w.j)) ', ' pctxt],'Interpreter','none','FontWeight','normal');
else
    title(ax,{[rowtxt ', ' char(ud.signames(w.j)) ', ' pctxt], 'grey = 1 s before/after, faded to full = the bout', 'dot = bout start, square = end, red star = stim'},'Interpreter','none','FontWeight','normal');
end
end
