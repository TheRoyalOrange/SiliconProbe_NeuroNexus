% ITPC_PhaseReset_singleanimal.m
%
% Description: For one animal and every probe in the recording (or the
%   ones listed in regions), pools peri-stimulus LFP trials from a
%   user-listed set of conditions into one reference condition
%   (ref_condition), then works on the LFP only (2026-11-02; CSD and the
%   earlier test plots are commented out in %{ %} blocks, kept for later):
%     1. parameters
%     2. loads ProbeInfo and picks the analysed channels on every probe
%        (probe channels not in that probe's Ch_Remove, with their probe,
%        row and shank in the probe map)
%     3. picks the ref_condition LFP files (uipickfiles), works out the
%        kept trials per file and loads their LFP segment for every
%        analysed channel
%     4 / 4b. OFF: per-trial inverse CSD via CSD_flex.m and its plot
%     5. 0.5-150 Hz bandpass (lfp_bp_hz) then complex Morlet wavelet on
%        lfp_ref (method from the OpenEphys_BaseAnalysis.m TF
%        section, cycles range_cycles, max(abs()) normalisation): phase per
%        channel x frequency x time x trial; power as trial means (raw and
%        per-trial-baseline dB). (The old CSD + LFP loop version is OFF.)
%     5b. mean pre vs post stim power spectra per channel, +-500 and
%        +-200 ms
%     6. per channel x frequency x time over ref_condition trials: ITPC,
%        preferred angle and Rayleigh Z; 2nd-harmonic ITPC and its axis
%        (catches trials split into two groups pi apart, where ITPC ~0);
%        and a descriptive unimodal/bimodal label map (phasemode_lfp).
%        (The old CSD + LFP loop version is OFF.)
%     7 / 8. OFF: CSD/LFP/difference test figures and the probe 1 vs
%        probe 2 volume-conduction check (PLV, lag, PSI)
%     9. per probe: pre (black) vs post (red) power spectra in the probe
%        layout, log power, one figure per window size
%     10. per probe: ITPC (itpcplot_ms, itpcplot_hz) in the probe layout
%        with Rayleigh outlines at p ~ 0.05 (cyan) and 0.01 (white),
%        thresholds from the trial count, uncorrected
%     10b. per probe: the same ITPC, blanked (NaN) outside Rayleigh
%        p ~ 0.05, so only the significant regions show, with the same
%        outlines and the stim-onset line
%     11. asks for a frequency band (inputdlg, repeats until you're happy),
%        bandpasses the trial-average LFP in that band and draws that ERP
%        over the section-10 ITPC plot on a right-hand y axis
%   Figures are not saved. Planned later: merge steps behind 0/1 toggles,
%   other conditions scored against ref_condition (e.g. a per-trial
%   phase-reset score), further connectivity measures.
%
% Inputs:
%   ProbeInfo (struct) - loaded from <probeinfo_dir>\<animal>-ProbeInfo.mat
%     .Areas (cell of char, 1 x probenum) - region per probe; regions picks
%       from these ({} = all)
%     .ProbeIds (cell, 1 x probenum) - raw channel IDs per probe (raw data
%       rows, NOT 1:N per probe - see oshkosh_detect_funct.m header)
%     .Ch_Remove (cell, 1 x probenum) - raw channel IDs to exclude
%     .ProbeMaps (cell, 1 x probenum) - 2D [depth rows x shanks] matrix of
%       raw channel IDs; columns = shanks. INFERRED: rows are equally
%       spaced contacts in physical order along the shank; which end row 1
%       is (superficial or deep) is not stated anywhere.
%   <animal>-<stim>-LFP.mat files from OpenEphys_BaseAnalysis.m,
%     chosen per condition with uipickfiles:
%     stim_lfp_stimchunks (double, trials x 15000 x channels) - LFP in
%       microvolts at 1 kHz. Epoch = stim-5000 : stim+9999 samples, so
%       stimulus onset is sample 5001 (onset_idx). Channel index = raw
%       channel ID.
%     stim_times (double, trials x 1) - stimulus onsets, 1 kHz samples
%     tr_keep (int vector) - kept trial indices
%     tr_remove - INFERRED, two formats exist: BaseAnalysis saves removed
%       trial indices (already excluded from tr_keep);
%       ConditionalTrialRemove_callable.m migrates missing/empty ones to a
%       0/1 mask the same length as tr_keep. Detected per file and logged
%       in report.files.tr_remove_format.
%     tr_remove_conditional (table, optional) - Name/TrialIdx rows from
%       ConditionalTrialRemove_callable.m; TrialIdx = original trial
%       numbers, NaN = condition checked, nothing removed. Applied only for
%       names listed in tr_conditional_names.
%   User parameters - see section 1.
%
% Outputs (workspace):
%   prbs (double, 1 x nprobes) - probe indices analysed (into the ProbeInfo
%     cell fields)
%   chansuse (double, nch x 1) - raw channel IDs analysed, every probe in
%     prbs stacked in order (each probe's channels minus its Ch_Remove).
%     Index k of every channel dimension = chansuse(k).
%   ch_probe, ch_row, ch_shank (double, nch x 1) - probe index, row (depth
%     position) and column (shank) of each chansuse channel in
%     ProbeInfo.ProbeMaps{ch_probe(k)}
%   lfp_ref (double, trials x nseg x nch) - ref_condition LFP in
%     microvolts, segment = time_ms padded by pad_ms on both sides
%   seg_ms (double, 1 x nseg) - time of each segment sample, ms re: stim
%   trial_table (table, one row per ref_condition trial, same order as
%     dim 1 of lfp_ref) - condition, condiname, file, filename, trial
%     (original trial number in its file), stim_time
%   phase_lfp (single, nch x nfrex x ntime x ntrials) - wavelet phase of
%     the bandpassed lfp_ref, radians (-pi, pi]. Trial order = trial_table
%     rows.
%   pow_lfp, powdB_lfp (single, nch x nfrex x ntime) - trial-mean wavelet
%     power (uV^2), and trial mean of per-trial dB vs tf_baseline_ms
%   spec_pre500, spec_post500, spec_pre200, spec_post200 (single, nch x
%     nfrex) - pow_lfp averaged over [-500 0), (0 500], [-200 0), (0 200]
%     ms, uV^2. At low frequencies the wavelet smears across 0 ms.
%   itpc_lfp, prefAngle_lfp, rayleighZ_lfp, itpc2_lfp, prefAxis_lfp,
%     phasemode_lfp - section 6, same meaning as the *_ref versions below
%     but from phase_lfp
%   erp_band (double, 1 x 2) - last band picked in section 11, Hz
%   erp_band_lfp (double, nch x nseg) - trial-average LFP bandpassed to
%     erp_band (IIR, zero-phase), uV, over seg_ms
%   --- the outputs below come from sections that are OFF for now ---
%   csd_ref (double, trials x nseg x nch) - ref_condition inverse CSD in
%     uA/mm^3, same dims/order as lfp_ref (positive = source, negative =
%     sink). NaN on shanks with fewer than csd_min_chans kept channels.
%   frex (double, 1 x nfrex) - wavelet frequencies, Hz
%   times_ms (double, 1 x ntime) - time of each TF sample, ms re: stim
%     (= time_ms(1):time_ms(2))
%   phase_ref (single, nch x nfrex x ntime x ntrials) - wavelet phase of
%     csd_ref, radians (-pi, pi]. Trial order = trial_table rows. NaN on
%     shanks with no CSD.
%   pow_ref (single, nch x nfrex x ntime) - trial-mean wavelet power,
%     (uA/mm^3)^2
%   powdB_ref (single, nch x nfrex x ntime) - trial mean of per-trial dB
%     power, 10*log10(power / mean power in tf_baseline_ms) per frequency
%   itpc_ref (single, nch x nfrex x ntime) - ITPC, |mean(exp(i*phase))|
%     over trials, 0-1. Biased upward for small n.
%   prefAngle_ref (single, nch x nfrex x ntime) - preferred phase,
%     angle(mean(exp(i*phase))), radians. Meaningless where ITPC ~ chance.
%   rayleighZ_ref (single, nch x nfrex x ntime) - ntrials * ITPC^2; ~1
%     under no phase locking whatever the trial count, so comparable
%     across conditions with different n
%   itpc2_ref (single, nch x nfrex x ntime) - 2nd-harmonic ITPC,
%     |mean(exp(2i*phase))|: high when trials split into two groups pi
%     apart. A tight unimodal cluster is also high here (wrapped normal:
%     ITPC2 = ITPC^4), so read it together with phasemode_ref.
%   prefAxis_ref (single, nch x nfrex x ntime) - axis of a pi-apart split,
%     angle(mean(exp(2i*phase)))/2, radians, defined modulo pi
%   phasemode_ref (int8, nch x nfrex x ntime) - 0 = no significant
%     clustering, 1 = unimodal (Rayleigh test on phase p < phasemode_alpha),
%     2 = bimodal (Rayleigh on 2*phase p < phasemode_alpha and unimodal
%     test not significant). Uncorrected, descriptive guide only: with no
%     phase locking ~2*phasemode_alpha of points are labelled by chance.
%   p12 (struct, section 8 only) - probe 1 (region) vs probe 2 (p2_region):
%     .plv, .lag, .lag_ms (single, nch1 x nch2 x nwin x 2 x 2; dims in
%       .dims = probe1 chan, probe2 chan, window, signal {LFP,CSD}, ERP
%       {kept,removed}). plv = mean over the window's freqs and times of
%       |trial-mean exp(i(phi1-phi2))|; lag = angle of the summed trial-mean
%       vectors, radians, positive = probe 1 phase ahead; lag_ms = lag at
%       the band centre frequency. NaN for channels on no-CSD shanks (CSD).
%     .psi, .psi_z (single, same size) - phase-slope index over the band's
%       1 Hz steps, from coherency (cross-spectrum over trials and window
%       times); psi_z = psi / jackknife SD over trials. + = probe 1 leads,
%       ~0 for volume conduction; |psi_z| > 2 ~ p < 0.05, uncorrected.
%     .plv_t, .lag_t (single, nch1 x nch2 x ntt x nwin x 2 x 2; dims in
%       .dims_t) - same PLV / lag (radians) but at every time point in
%       p2_tr_ms, averaged over the band only, so brief coupling isn't
%       diluted by the window average. .times_tr (1 x ntt) - ms re: stim.
%     .chans1 (= chansuse), .chans2 (probe 2 raw IDs minus its Ch_Remove),
%     .region1, .region2, .windows (p2_windows), .ntrials
%   report (struct) - run bookkeeping for the end-of-run report: channel
%     counts and removed channels per probe (indexed by probe), per-file
%     trial table (report.files), phase_lfp size (report.phase_gb), notes
%
% Dependencies: OpenEphys_BaseAnalysis.m (LFP files, ProbeInfo),
%   OpenEphys_editProbeInfo_ChRemove.m (Ch_Remove), uipickfiles.m,
%   Signal Processing Toolbox (bandpass); optional
%   ConditionalTrialRemove_callable.m (tr_remove_conditional).
%   OFF sections also need CSD_flex.m (this folder) and slanCM.m.

%% 1 PARAMETERS
animal = '20261011-p7';
%region = 'V1'; %must match one entry of ProbeInfo.Areas. Old: one probe per run, now every probe in one go (regions below)

%which probes do you want? names from ProbeInfo.Areas, e.g. {'V1','S1'}
regions = {}; %{} = every probe in the recording

ref_condis = {'whiskerActive_TTX'}; %conditions pooled into the reference condition (one file-picker prompt each)
ref_condition = 'WhiskerActive_TTX'; %name of the pooled reference condition (titles, outputs)

time_ms = [-1000 1000]; %analysis window, ms re: stim
%pad_ms = 1000; %old: enough for 2 Hz, not for the 1 Hz wavelet + 0.5 Hz highpass now in use
pad_ms = 1500; %extra data loaded on each side so wavelet and filter edge effects stay outside time_ms (1 Hz wavelet reaches ~1.4 s)
onset_idx = 5001; %stimulus onset sample in stim_lfp_stimchunks

tr_conditional_names = {}; %e.g. {'whisker_twitch_artifact'}; extra trial exclusions

%inverse CSD (CSD_flex.m; defaults = the CSDallshank call in OpenEphys_BaseAnalysis.m)
%CSD is off for now (LFP only), so these are parked too
%csd_spacing_um = 100; %distance between adjacent ProbeMaps rows, um
%csd_radius_um = 300; %lateral radius of the active column of tissue, um (not the contact size: contacts are 177 um^2 in area, NeuroNexus catalog A8x8-Edge-5mm-100-200-177)
%csd_conductivity = 0.3; %extracellular conductivity, S/m
%csd_min_chans = 3; %shanks with fewer kept channels get NaN CSD

%broadband clean-up before the wavelet (the LFP files are only downsampled, never filtered)
lfp_bp_hz = [0.5 150]; %Hz, bandpass edges

%wavelet (method from the OpenEphys_BaseAnalysis.m TF section)
%frex = 2:1:80; %old range, now down to 1 Hz for the 1-80 Hz ITPC plot
frex = 1:1:80; %Hz, integer steps so frequency index is unambiguous
range_cycles = [3 5]; %wavelet cycles at lowest/highest frequency (BaseAnalysis uses [6 10]; fewer cycles = better time resolution)
tf_baseline_ms = [-800 -200]; %baseline window for dB power, ms re: stim (must be inside time_ms)
max_phase_gb = 30; %stop before allocating phase_lfp if it would be larger than this (GB)

%clear the per-trial phases (phase_lfp) once ITPC is done? it's the biggest thing in memory (~16 GB for 128 channels x 200 trials)
mr_clean = 0;
    % 0 = keep phase_lfp (e.g. for per-trial phase-reset scores later)
    % 1 = clear phase_lfp at the end of section 6

%phase clustering (section 6)
phasemode_alpha = 0.01; %uncorrected Rayleigh threshold for the phasemode_ref label map (descriptive, not inference)

%%%Plotting details%%%

%what time and frequency range for the ITPC plots?
itpcplot_ms = [-200 500]; %ms re: stim
itpcplot_hz = [1 80]; %Hz
itpc_clim = [0 0.4]; %ITPC colour limits

%which band should the band-ERP section offer first? (you get asked anyway)
erp_band = [4 8]; %Hz

probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';

%% 2 LOAD PROBEINFO 
report = struct();
report.notes = {};

load(fullfile(probeinfo_dir, [animal '-ProbeInfo.mat']), 'ProbeInfo');

%old: one probe (region) per run. Replaced by the every-probe version below
%{
prb = find(ismember(ProbeInfo.Areas, region));
if numel(prb) ~= 1
    error('Region %s matches %d probes in ProbeInfo.Areas', region, numel(prb));
end

%channels: raw IDs on this probe minus Ch_Remove. Index k of every channel dim = chansuse(k)
chansall = ProbeInfo.ProbeIds{prb}(:);
chremove = [];
if numel(ProbeInfo.Ch_Remove) >= prb
    chremove = ProbeInfo.Ch_Remove{prb};
end
chansuse = chansall(~ismember(chansall, chremove));
nch = numel(chansuse);

probemap = ProbeInfo.ProbeMaps{prb}; %[depth rows x shanks]
ch_row = zeros(nch,1);
ch_shank = zeros(nch,1);
for ch = 1:nch
    [ch_row(ch), ch_shank(ch)] = find(probemap == chansuse(ch));
end
report.nchans_probe = numel(chansall);
report.nchans_used = nch;
report.ch_removed = chremove(:)';
%}

%which probes are we doing?
if isempty(regions) == true
    prbs = 1:numel(ProbeInfo.Areas); %nothing picked = every probe
else
    prbs = find(ismember(ProbeInfo.Areas, regions)); %just the ones you named
end
if isempty(prbs) == true
    error('None of regions {%s} are in ProbeInfo.Areas {%s}', strjoin(regions, ', '), strjoin(ProbeInfo.Areas, ', '));
end

%channels: raw IDs on each probe minus its Ch_Remove, stacked probe after probe.
%Index k of every channel dim = chansuse(k), and ch_probe(k) says which probe it's on
chansuse = [];
ch_probe = [];
ch_row = [];
ch_shank = [];
for prb = prbs
    chansall = ProbeInfo.ProbeIds{prb}(:); %raw channel IDs on this probe
    chremove = [];
    if numel(ProbeInfo.Ch_Remove) >= prb
        chremove = ProbeInfo.Ch_Remove{prb}; %channels marked bad on this probe
    end
    chans_p = chansall(~ismember(chansall, chremove)); %the keepers

    probemap = ProbeInfo.ProbeMaps{prb}; % [depth rows x shanks]
    row_p = zeros(numel(chans_p),1);
    shank_p = zeros(numel(chans_p),1);
    for ch = 1:numel(chans_p)
        [row_p(ch), shank_p(ch)] = find(probemap == chans_p(ch)); %where this channel sits on the probe
    end

    chansuse = [chansuse; chans_p]; %#ok<AGROW>
    ch_probe = [ch_probe; repmat(prb, numel(chans_p), 1)]; %#ok<AGROW>
    ch_row = [ch_row; row_p]; %#ok<AGROW>
    ch_shank = [ch_shank; shank_p]; %#ok<AGROW>
    report.nchans_probe(prb) = numel(chansall); %bookkeeping for the end-of-run report
    report.nchans_used(prb) = numel(chans_p);
    report.ch_removed{prb} = chremove(:)';
end
nch = numel(chansuse); %all probes together


%% 3 LFP FILES FOR REF_CONDITION
%pick LFP files per condition
condi_dir = cell(1,numel(ref_condis));
for con = 1:numel(ref_condis)
    pickedfiles = uipickfiles('FilterSpec', fullfile(lfp_dir, animal), 'Prompt', ['Choose ' ref_condis{con} ' LFP files (pooled into ' ref_condition ')']);
    if ~iscell(pickedfiles) || isempty(pickedfiles)
        error('No files chosen for condition %s', ref_condis{con});
    end
    condi_dir{con} = pickedfiles;
end

%segment loaded per trial
seg_ms = (time_ms(1)-pad_ms):(time_ms(2)+pad_ms);
seg_idx = onset_idx + seg_ms;
nseg = numel(seg_idx);
if seg_idx(1) < 1
    error('time_ms(1) - pad_ms reaches before the start of the epoch.');
end

%trials used per file (tr_keep, minus tr_remove in either format, minus chosen conditional removals), and their LFP
data_lfp = {};
tt = {};
filerows = {};
for con = 1:numel(ref_condis)
    for file = 1:numel(condi_dir{con})
        fpath = condi_dir{con}{file};
        [~, fname] = fileparts(fpath);
        vars = who('-file', fpath);
        F = load(fpath, 'tr_keep', 'tr_remove', 'stim_times');
        tr_keep = F.tr_keep(:)';
        tr_remove = [];
        if isfield(F, 'tr_remove')
            tr_remove = F.tr_remove(:)';
        end
        if ~isempty(tr_remove) && numel(tr_remove) == numel(tr_keep) && all(ismember(tr_remove, [0 1]))
            tr_use = tr_keep(~logical(tr_remove)); %0/1 mask over tr_keep (ConditionalTrialRemove format)
            rmformat = "mask";
        else
            tr_use = setdiff(tr_keep, tr_remove, 'stable'); %index list (BaseAnalysis format)
            rmformat = "index";
        end

        ncond_removed = 0;
        if ~isempty(tr_conditional_names)
            if ismember('tr_remove_conditional', vars)
                C = load(fpath, 'tr_remove_conditional');
                trc = C.tr_remove_conditional;
                missing = setdiff(string(tr_conditional_names), string(trc.Name));
                if ~isempty(missing)
                    report.notes{end+1} = sprintf('%s: conditional name(s) not recorded: %s', fname, strjoin(missing, ', '));
                end
                rm = trc.TrialIdx(ismember(string(trc.Name), string(tr_conditional_names)));
                rm = rm(~isnan(rm)); %NaN = checked, nothing removed
                ncond_removed = sum(ismember(tr_use, rm));
                tr_use = tr_use(~ismember(tr_use, rm));
            else
                report.notes{end+1} = sprintf('%s: no tr_remove_conditional table, conditional removals skipped', fname);
            end
        end

        %LFP for the kept trials, segment and channels (file is not -v7.3, so load it whole once)
        S = load(fpath, 'stim_lfp_stimchunks');
        sz = size(S.stim_lfp_stimchunks);
        if sz(2) < seg_idx(end) || sz(3) < max(chansuse) || sz(1) < max([tr_use 0])
            error('%s: stim_lfp_stimchunks is %s, too small for this segment/channel/trial set', fname, mat2str(sz));
        end
        data_lfp{end+1} = double(S.stim_lfp_stimchunks(tr_use, seg_idx, chansuse)); %#ok<SAGROW> %trials x nseg x nch
        clear S

        n = numel(tr_use);
        tt{end+1} = table(repmat(con,n,1), repmat(string(ref_condis{con}),n,1), repmat(file,n,1), ...
            repmat(string(fname),n,1), tr_use(:), F.stim_times(tr_use(:)), ...
            'VariableNames', {'condition','condiname','file','filename','trial','stim_time'}); %#ok<SAGROW>
        filerows(end+1,:) = {string(ref_condis{con}), file, string(fname), numel(F.stim_times), numel(tr_keep), rmformat, ncond_removed, n}; %#ok<SAGROW>
    end
end
lfp_ref = cat(1, data_lfp{:}); %trials x nseg x nch, microvolts
clear data_lfp
trial_table = vertcat(tt{:});
report.files = cell2table(filerows, 'VariableNames', ...
    {'condition','file','filename','ntrials_total','ntrials_keep','tr_remove_format','ncond_removed','ntrials_used'});
report.ntrials = size(lfp_ref,1);

%% 4 Data PreProcessing: per-trial inverse CSD along each shank (OFF FOR NOW)
%CSD switched off 2026-11-02: sticking to the LFP for now. Code kept for later, just take away the %{ %} around it
%{
%lfp_ref holds chansuse (Ch_Remove already dropped), so 'ChRemove' is true
[csd_ref, report.csd] = CSD_flex(lfp_ref, ProbeInfo, prb, 'ChRemove', true, ...
    'spacing', csd_spacing_um, 'radius', csd_radius_um, 'conductivity', csd_conductivity, 'minchans', csd_min_chans);
if ~isequal(report.csd.chan_ids, chansuse)
    error('CSD_flex channel order does not match chansuse.');
end

%}

%% 4b Data Visualization: CSD plot (OFF FOR NOW)
%off along with the rest of the CSD stuff
%{
%trial-mean CSD per shank, layout borrowed from the OpenEphys_BaseAnalysis.m CSD plot.
%Removed (Ch_Remove) rows stay empty at their true depth (dashed grey line).
csdplot_ms = [-50 500]; %plot window, ms re: stim
csdplot_overlay = 'LFP'; %'LFP', 'CSD' or 'none': normalised trial-mean trace drawn on each row
csdplot_cmap = slanCM('guppy');

pidx = find(seg_ms >= csdplot_ms(1) & seg_ms <= csdplot_ms(2));
ptime = seg_ms(pidx);
csd_avg = reshape(mean(csd_ref(:,pidx,:), 1, 'omitnan'), numel(pidx), []); %time x nch
lfp_avg = reshape(mean(lfp_ref(:,pidx,:), 1), numel(pidx), []); %time x nch
csdmap = ProbeInfo.ProbeMaps{prb}; %[depth rows x shanks]
[nrow, nshank] = size(csdmap);
lim = 1.1*max(abs(csd_avg(:)), [], 'omitnan'); %one symmetric scale shared by all shanks
if isempty(lim) || isnan(lim) || lim == 0
    lim = 1;
end

figure('Position', get(0, 'Screensize'));
sgtitle([animal ' ' region '-Probe ' ref_condition ' CSD (mean of ' num2str(size(csd_ref,1)) ' trials)'], 'Interpreter', 'none')
for q = 1:nshank
    [iskept, loc] = ismember(csdmap(:,q), chansuse);
    shank_csd = nan(nrow, numel(pidx)); %rows x time, NaN for removed rows
    shank_csd(iskept,:) = csd_avg(:, loc(iskept))';

    subplot(1, nshank, q)
    if ~all(isnan(shank_csd(:))) %contour warns on an all-NaN shank
        contour(ptime, 1:nrow, shank_csd, 200);
    end
    set(gca, 'Color', [0 0 0], 'YDir', 'reverse', 'CLim', [-lim lim])
    ylim([.5 nrow+.5])
    xlim([ptime(1) ptime(end)])
    colormap(gca, csdplot_cmap)
    c = colorbar('southoutside');
    c.Label.String = 'CSD (uA/mm^3)';
    if all(isnan(shank_csd(:)))
        title(['Shank ' num2str(q) ' (no CSD: < csd_min_chans)'], 'Interpreter', 'none')
    else
        title(['Shank ' num2str(q)])
    end
    xlabel('Time (ms)')
    ylabel('ProbeMaps row')

    hold on
    for j = 1:nrow
        if ~iskept(j)
            yline(j, '--', 'Color', [.5 .5 .5], 'LineWidth', 1);
            continue
        end
        switch csdplot_overlay
            case 'LFP'
                tr = lfp_avg(:, loc(j));
            case 'CSD'
                tr = csd_avg(:, loc(j));
            otherwise
                tr = [];
        end
        if ~isempty(tr) && any(tr ~= 0) && ~any(isnan(tr))
            plot(ptime, -tr/max(abs(tr))/2 + j, 'Color', 'white', 'LineWidth', .5);
        end
    end
    xline(0, 'Color', 'g', 'LineWidth', 2);
    hold off
end

%}

%% 5 (OFF) old wavelet: CSD + LFP loop
%replaced 2026-11-02 by the LFP-only section 5 just below (no CSD/LFP loop). Kept for when the CSD comes back
%{
%same wavelet, trials and time points for both signals, so CSD- and LFP-based results compare directly
times_ms = time_ms(1):time_ms(2);
ntime = numel(times_ms);
nfrex = numel(frex);
ntr = size(csd_ref,1);
[~, keep_idx] = ismember(times_ms, seg_ms); %positions of time_ms inside the padded segment
bidx = find(times_ms >= tf_baseline_ms(1) & times_ms <= tf_baseline_ms(2)); %baseline positions in times_ms
if tf_baseline_ms(1) < time_ms(1) || tf_baseline_ms(2) > time_ms(2) || isempty(bidx)
    error('tf_baseline_ms [%d %d] must lie inside time_ms [%d %d].', tf_baseline_ms, time_ms);
end
report.phase_gb = 2*4*nch*nfrex*ntime*ntr/1e9; %phase_ref + phase_lfp
if report.phase_gb > max_phase_gb
    error('phase_ref + phase_lfp would need %.1f GB (> max_phase_gb = %g). Narrow frex or time_ms, or raise max_phase_gb.', report.phase_gb, max_phase_gb);
end

%wavelet parameters (as in OpenEphys_BaseAnalysis.m TF section)
s = logspace(log10(range_cycles(1)), log10(range_cycles(end)), nfrex) ./ (2*pi*frex); %Gaussian width, s
wavtime = -2:1/1000:2;
half_wave = (length(wavtime)-1)/2;
nWave = length(wavtime);
if 3*s(1)*1000 > pad_ms
    report.notes{end+1} = sprintf(['Lowest-frequency wavelet (%g Hz, 3 SD = %.0f ms) reaches past pad_ms = %d ms, ' ...
        'so edge effects reach into time_ms at that frequency.'], frex(1), 3*s(1)*1000, pad_ms);
end

%wavelet FFTs, built once: all trials of a channel go through one convolution
nConv = nWave + nseg*ntr - 1;
waveletX = complex(zeros(nfrex, nConv, 'single'));
for fi = 1:nfrex
    wavelet = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    wX = fft(wavelet, nConv);
    %waveletX(fi,:) = single(wX ./ max(wX)); %BaseAnalysis version: max(wX) is complex, so it also rotates every phase by a constant that depends on nConv (i.e. on trial count)
    waveletX(fi,:) = single(wX ./ max(abs(wX))); %normalise gain only, phase untouched
end

for src = 1:2 %1 = CSD (csd_ref -> *_ref outputs), 2 = LFP (lfp_ref -> *_lfp outputs)
    ph_src = nan(nch, nfrex, ntime, ntr, 'single');
    pow_src = nan(nch, nfrex, ntime, 'single');
    powdB_src = nan(nch, nfrex, ntime, 'single');
    for ch = 1:nch
        if src == 1
            dat = csd_ref(:,:,ch); %trials x nseg
        else
            dat = lfp_ref(:,:,ch);
        end
        if all(isnan(dat(:)))
            continue %shank with no CSD (< csd_min_chans)
        end
        %demean each trial: the LFP is unfiltered (BaseAnalysis only downsamples) and the CSD keeps
        %per-channel offsets; 3-cycle wavelets pass ~1% of a DC offset, which would add the same
        %phase to every trial and inflate low-frequency ITPC
        dat = dat - mean(dat, 2);
        dataX = fft(reshape(dat', 1, []), nConv); %trials concatenated end to end

        for fi = 1:nfrex
            as = ifft(double(waveletX(fi,:)) .* dataX);
            as = as(half_wave+1:end-half_wave);
            as = reshape(as, nseg, ntr); %time x trials
            as = as(keep_idx,:);
            ph_src(ch,fi,:,:) = reshape(single(angle(as)), [1 1 ntime ntr]);
            pw = abs(as).^2;
            pow_src(ch,fi,:) = single(mean(pw, 2));
            powdB_src(ch,fi,:) = single(mean(10*log10(pw ./ mean(pw(bidx,:), 1)), 2));
        end
    end
    if src == 1
        phase_ref = ph_src; pow_ref = pow_src; powdB_ref = powdB_src;
    else
        phase_lfp = ph_src; pow_lfp = pow_src; powdB_lfp = powdB_src;
    end
    clear ph_src pow_src powdB_src
end
clear waveletX wavelet wX dataX as pw dat

%}



%Now beginning the time-frequency section (LFP only)

%% 5 Data PreProcessing: bandpass + Morlet wavelet on the LFP (phase per trial, power)
%same wavelet method as the OpenEphys_BaseAnalysis.m TF section (logspace cycles, wavtime -2:2,
%all trials of a channel convolved in one go), with range_cycles [3 5] and a padded segment
disp('TF: Bandpassing and wavelet-ing the LFP, sit tight')
tic
times_ms = time_ms(1):time_ms(2); %ms re: stim, the part we keep after convolution
ntime = numel(times_ms);
nfrex = numel(frex);
ntr = size(lfp_ref,1); %number of trials
[~, keep_idx] = ismember(times_ms, seg_ms); %where time_ms sits inside the padded segment
bidx = find(times_ms >= tf_baseline_ms(1) & times_ms <= tf_baseline_ms(2)); %baseline points for dB power
if tf_baseline_ms(1) < time_ms(1) || tf_baseline_ms(2) > time_ms(2) || isempty(bidx)
    error('tf_baseline_ms [%d %d] must lie inside time_ms [%d %d].', tf_baseline_ms, time_ms);
end
report.phase_gb = 4*nch*nfrex*ntime*ntr/1e9; %how big phase_lfp will be (single = 4 bytes)
if report.phase_gb > max_phase_gb
    error('phase_lfp would need %.1f GB (> max_phase_gb = %g). Narrow frex or time_ms, or raise max_phase_gb.', report.phase_gb, max_phase_gb);
end

%wavelet parameters (as in OpenEphys_BaseAnalysis.m TF section)
s = logspace(log10(range_cycles(1)), log10(range_cycles(end)), nfrex) ./ (2*pi*frex); %Gaussian width, s
wavtime = -2:1/1000:2;
half_wave = (length(wavtime)-1)/2;
nWave = length(wavtime);
if 3*s(1)*1000 > pad_ms
    report.notes{end+1} = sprintf(['Lowest-frequency wavelet (%g Hz, 3 SD = %.0f ms) reaches past pad_ms = %d ms, ' ...
        'so edge effects reach into time_ms at that frequency.'], frex(1), 3*s(1)*1000, pad_ms);
end

%wavelet FFTs, built once: all trials of a channel go through one convolution
nConv = nWave + nseg*ntr - 1;
waveletX = complex(zeros(nfrex, nConv, 'single')); % [frequencies x convolution length]
for fi = 1:nfrex
    wavelet = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2)); %create wavelet
    wX = fft(wavelet, nConv); %and get its FFT
    %waveletX(fi,:) = single(wX ./ max(wX)); %BaseAnalysis version: max(wX) is complex, so it also rotates every phase by a constant that depends on nConv (i.e. on trial count)
    waveletX(fi,:) = single(wX ./ max(abs(wX))); %normalise gain only, phase untouched
end

phase_lfp = zeros(nch, nfrex, ntime, ntr, 'single'); % [channels x frequencies x time x trials], radians
pow_lfp = zeros(nch, nfrex, ntime, 'single'); % [channels x frequencies x time], trial-mean power, uV^2
powdB_lfp = zeros(nch, nfrex, ntime, 'single'); % [channels x frequencies x time], trial mean of dB vs baseline
%for ch = 1:nch %old: one channel after another (~3-4 min for 128 channels x 200 trials)
parfor ch = 1:nch %channels don't depend on each other, so run them side by side (each worker gets its own copy of waveletX, ~0.6 GB)
    dat = lfp_ref(:,:,ch); % [trials x segment], uV
    %dat = dat - mean(dat, 2); %old: per-trial demean. The broadband bandpass below takes out DC and drift instead
    dat = bandpass(dat', lfp_bp_hz, 1000, 'ImpulseResponse', 'iir')'; %0.5-150 Hz clean-up, each trial filtered on its own (bandpass works down columns)
    dataX = fft(reshape(dat', 1, []), nConv); %trials concatenated end to end

    %this channel's results first, then they go into the big arrays in one go (parfor needs it that way)
    ph_ch = zeros(nfrex, ntime, ntr, 'single'); % [frequencies x time x trials]
    pow_ch = zeros(nfrex, ntime, 'single'); % [frequencies x time]
    powdB_ch = zeros(nfrex, ntime, 'single');
    for fi = 1:nfrex
        as = ifft(double(waveletX(fi,:)) .* dataX); %run convolution in one step
        as = as(half_wave+1:end-half_wave); %trim off the wavelet's overhang
        as = reshape(as, nseg, ntr); % [time x trials]
        as = as(keep_idx,:); %just time_ms
        ph_ch(fi,:,:) = reshape(single(angle(as)), [1 ntime ntr]); %phase of every trial
        pw = abs(as).^2; %power
        pow_ch(fi,:) = single(mean(pw, 2)); %trial-mean power
        powdB_ch(fi,:) = single(mean(10*log10(pw ./ mean(pw(bidx,:), 1)), 2)); %each trial vs its own baseline, then averaged
    end
    phase_lfp(ch,:,:,:) = reshape(ph_ch, [1 nfrex ntime ntr]);
    pow_lfp(ch,:,:) = reshape(pow_ch, [1 nfrex ntime]);
    powdB_lfp(ch,:,:) = reshape(powdB_ch, [1 nfrex ntime]);
    disp(['TF: channel ' num2str(ch) ' of ' num2str(nch) ' done (workers finish in any order)'])
end
clear waveletX wavelet wX dataX as pw dat
elapsed_time = toc;
fprintf('execution time: %.2f seconds\n', elapsed_time);

%% 5b Data PreProcessing: pre vs post stim power spectra (mean power in each window)
%heads up: at low frequencies the wavelet is long (at 1 Hz it smears ~1 s each way), so pre and post bleed into each other there
pre500 = times_ms >= -500 & times_ms < 0; %which time points are in each window
post500 = times_ms > 0 & times_ms <= 500;
pre200 = times_ms >= -200 & times_ms < 0;
post200 = times_ms > 0 & times_ms <= 200;
pre100 = times_ms >= -100 & times_ms < 0;
post100 = times_ms > 0 & times_ms <= 100;
spec_pre500 = mean(pow_lfp(:,:,pre500), 3); % [channels x frequencies], uV^2
spec_post500 = mean(pow_lfp(:,:,post500), 3);
spec_pre200 = mean(pow_lfp(:,:,pre200), 3);
spec_post200 = mean(pow_lfp(:,:,post200), 3);
spec_pre100 = mean(pow_lfp(:,:,pre100), 3);
spec_post100 = mean(pow_lfp(:,:,post100), 3);

%% 6 (OFF) old ITPC: CSD + LFP loop
%replaced 2026-11-02 by the LFP-only section 6 just below (no CSD/LFP loop). Kept for when the CSD comes back
%{
%one channel at a time: exp(1i*phase) on the whole array would be a complex double copy 4x its size
ntr = size(phase_ref,4);
nfrex = size(phase_ref,2);
ntime = size(phase_ref,3);

%Rayleigh test p-value (Zar 1999, accurate for small n); Rn = n * resultant length
rayleigh_p = @(R) min(1, exp(sqrt(1 + 4*ntr + 4*(ntr^2 - (ntr*R).^2)) - (1 + 2*ntr)));

for src = 1:2 %1 = CSD (phase_ref -> *_ref), 2 = LFP (phase_lfp -> *_lfp)
    itpc_src = nan(nch, nfrex, ntime, 'single');
    prefAngle_src = nan(nch, nfrex, ntime, 'single');
    rayleighZ_src = nan(nch, nfrex, ntime, 'single');
    itpc2_src = nan(nch, nfrex, ntime, 'single');
    prefAxis_src = nan(nch, nfrex, ntime, 'single');
    phasemode_src = zeros(nch, nfrex, ntime, 'int8');
    for ch = 1:nch
        if src == 1
            ph = double(reshape(phase_ref(ch,:,:,:), nfrex, ntime, ntr));
        else
            ph = double(reshape(phase_lfp(ch,:,:,:), nfrex, ntime, ntr));
        end
        if all(isnan(ph(:)))
            continue %shank with no CSD
        end
        z1 = mean(exp(1i*ph), 3); %first harmonic: one preferred angle
        z2 = mean(exp(2i*ph), 3); %second harmonic: two groups pi apart fold onto each other
        R1 = abs(z1);
        R2 = abs(z2);
        itpc_src(ch,:,:) = R1;
        prefAngle_src(ch,:,:) = angle(z1);
        rayleighZ_src(ch,:,:) = ntr*R1.^2;
        itpc2_src(ch,:,:) = R2;
        prefAxis_src(ch,:,:) = angle(z2)/2;

        %label: unimodal if phase clusters; bimodal only if 2*phase clusters and phase does not
        %(a tight unimodal cluster also has high ITPC2)
        uni = rayleigh_p(R1) < phasemode_alpha;
        bi = rayleigh_p(R2) < phasemode_alpha & ~uni;
        mode_ch = zeros(nfrex, ntime, 'int8');
        mode_ch(uni) = 1;
        mode_ch(bi) = 2;
        phasemode_src(ch,:,:) = mode_ch;
    end
    if src == 1
        itpc_ref = itpc_src; prefAngle_ref = prefAngle_src; rayleighZ_ref = rayleighZ_src;
        itpc2_ref = itpc2_src; prefAxis_ref = prefAxis_src; phasemode_ref = phasemode_src;
    else
        itpc_lfp = itpc_src; prefAngle_lfp = prefAngle_src; rayleighZ_lfp = rayleighZ_src;
        itpc2_lfp = itpc2_src; prefAxis_lfp = prefAxis_src; phasemode_lfp = phasemode_src;
    end
end
clear ph z1 z2 R1 R2 uni bi mode_ch itpc_src prefAngle_src rayleighZ_src itpc2_src prefAxis_src phasemode_src

%}

%% 6 Data PreProcessing: ITPC, preferred angle, Rayleigh stats and phase mode (LFP)
%one channel at a time: exp(1i*phase) on the whole array would be a complex double copy 4x its size
ntr = size(phase_lfp,4);
nfrex = size(phase_lfp,2);
ntime = size(phase_lfp,3);

%Rayleigh test p-value (Zar 1999, accurate for small n); Rn = n * resultant length
rayleigh_p = @(R) min(1, exp(sqrt(1 + 4*ntr + 4*(ntr^2 - (ntr*R).^2)) - (1 + 2*ntr)));

itpc_lfp = zeros(nch, nfrex, ntime, 'single'); % [channels x frequencies x time], 0-1
prefAngle_lfp = zeros(nch, nfrex, ntime, 'single'); %radians
rayleighZ_lfp = zeros(nch, nfrex, ntime, 'single'); %n * ITPC^2
itpc2_lfp = zeros(nch, nfrex, ntime, 'single'); %2nd harmonic ITPC
prefAxis_lfp = zeros(nch, nfrex, ntime, 'single'); %radians, modulo pi
phasemode_lfp = zeros(nch, nfrex, ntime, 'int8'); %0 none, 1 unimodal, 2 bimodal
for ch = 1:nch
    ph = double(reshape(phase_lfp(ch,:,:,:), nfrex, ntime, ntr)); % [frequencies x time x trials]
    z1 = mean(exp(1i*ph), 3); %first harmonic: one preferred angle
    z2 = mean(exp(2i*ph), 3); %second harmonic: two groups pi apart fold onto each other
    R1 = abs(z1);
    R2 = abs(z2);
    itpc_lfp(ch,:,:) = R1;
    prefAngle_lfp(ch,:,:) = angle(z1);
    rayleighZ_lfp(ch,:,:) = ntr*R1.^2;
    itpc2_lfp(ch,:,:) = R2;
    prefAxis_lfp(ch,:,:) = angle(z2)/2;

    %label: unimodal if phase clusters; bimodal only if 2*phase clusters and phase does not
    %(a tight unimodal cluster also has high ITPC2)
    uni = rayleigh_p(R1) < phasemode_alpha;
    bi = rayleigh_p(R2) < phasemode_alpha & ~uni;
    mode_ch = zeros(nfrex, ntime, 'int8');
    mode_ch(uni) = 1;
    mode_ch(bi) = 2;
    phasemode_lfp(ch,:,:) = mode_ch;
end
clear ph z1 z2 R1 R2 uni bi mode_ch

%done with the per-trial phases? they're the biggest thing in memory, and nothing after this section uses them
if mr_clean == 1 %set in section 1
    clear phase_lfp
end

%% 7 (OFF) test plots: CSD / LFP / CSD-LFP time-frequency figures
%switched off 2026-11-02: these were for testing; sections 9-11 below are the plots for now. Kept for later
%{
%one figure per map, one panel per channel at its ProbeMaps (row, shank) position.
%Style from the OpenEphys_BaseAnalysis.m TF plot. The older scripts build the grid
%from the transposed map, which only matches square (8x8) maps; this uses rows x shanks directly.
tfplot_ms = [-50 200]; %time window, ms re: stim
tfplot_hz = [frex(1) frex(end)]; %frequency range, Hz
itpc_clim = [0 0.4]; %ITPC colour limits (as in the old Plotting_ITPC_probe_singleanimal.m, removed 2026-11-07)
powdB_clim = []; %dB colour limits; [] = symmetric +/- max |dB| over all plotted panels (CSD and LFP)
rayZ_clim = []; %Rayleigh Z colour limits; [] = [0, max Z over all plotted panels (CSD and LFP)]
itpcdiff_clim = []; %CSD - LFP ITPC colour limits; [] = symmetric +/- max |difference|
rayZdiff_clim = []; %CSD - LFP Rayleigh Z colour limits; [] = symmetric +/- max |difference|
tfplot_levels = 50; %contour levels

tpi = find(times_ms >= tfplot_ms(1) & times_ms <= tfplot_ms(2));
fpi = find(frex >= tfplot_hz(1) & frex <= tfplot_hz(2));
tfmap = ProbeInfo.ProbeMaps{prb}; %[depth rows x shanks]
[nrow, nshank] = size(tfmap);
ntr_plot = height(trial_table);
itpc_diff = itpc_ref - itpc_lfp; %positive = stronger phase locking in the CSD
rayZ_diff = rayleighZ_ref - rayleighZ_lfp;

%auto colour limits, shared by the CSD and LFP versions of each figure
pick = @(A) reshape(A(:,fpi,tpi), [], 1);
if isempty(powdB_clim)
    powdB_clim = [-1 1]*max(abs([pick(powdB_ref); pick(powdB_lfp)]), [], 'omitnan');
end
if isempty(rayZ_clim)
    rayZ_clim = [0 max([pick(rayleighZ_ref); pick(rayleighZ_lfp)], [], 'omitnan')];
end
if isempty(itpcdiff_clim)
    itpcdiff_clim = [-1 1]*max(abs(pick(itpc_diff)), [], 'omitnan');
end
if isempty(rayZdiff_clim)
    rayZdiff_clim = [-1 1]*max(abs(pick(rayZ_diff)), [], 'omitnan');
end
cl = {powdB_clim, rayZ_clim, itpcdiff_clim, rayZdiff_clim};
for i = 1:numel(cl)
    if any(isnan(cl{i})) || diff(cl{i}) <= 0
        cl{i} = [0 1]; %all-NaN or all-zero data
    end
end
[powdB_clim, rayZ_clim, itpcdiff_clim, rayZdiff_clim] = cl{:};

%Rayleigh Z at p = 0.05 and 0.01 for this trial count: the Zar formula used in section 6,
%p = exp(sqrt(1 + 4n + 4(n^2 - Rn^2)) - (1 + 2n)), solved for Z = Rn^2/n (per point, uncorrected)
rayZ_p = [0.05 0.01];
L = log(rayZ_p) + 1 + 2*ntr_plot;
rayZ_thresh = (1 + 4*ntr_plot + 4*ntr_plot^2 - L.^2) / (4*ntr_plot);
rayZ_colors = {'c', 'w'}; %outline colour per rayZ_p level
masktitle = sprintf('ITPC where Rayleigh p < 0.05 (cyan p~0.05: Z=%.2f, white p~0.01: Z=%.2f, uncorrected)', rayZ_thresh);

%map, title, colormap, colour limits, colorbar label, Rayleigh Z map for outlines + mask ([] = none)
plotmaps = {itpc_ref, 'CSD ITPC', parula(256), itpc_clim, 'ITPC', []; ...
    powdB_ref, 'CSD Power (dB re: baseline)', slanCM('iceburn'), powdB_clim, 'dB', []; ...
    rayleighZ_ref, 'CSD Rayleigh Z', hot(256), rayZ_clim, 'Rayleigh Z', []; ...
    itpc_ref, ['CSD ' masktitle], parula(256), itpc_clim, 'ITPC', rayleighZ_ref; ...
    itpc_lfp, 'LFP ITPC', parula(256), itpc_clim, 'ITPC', []; ...
    powdB_lfp, 'LFP Power (dB re: baseline)', slanCM('iceburn'), powdB_clim, 'dB', []; ...
    rayleighZ_lfp, 'LFP Rayleigh Z', hot(256), rayZ_clim, 'Rayleigh Z', []; ...
    itpc_lfp, ['LFP ' masktitle], parula(256), itpc_clim, 'ITPC', rayleighZ_lfp; ...
    itpc_diff, 'ITPC CSD - LFP (positive = stronger in CSD)', slanCM('coolwarm'), itpcdiff_clim, 'ITPC difference', []; ...
    rayZ_diff, 'Rayleigh Z CSD - LFP (positive = stronger in CSD)', slanCM('coolwarm'), rayZdiff_clim, 'Rayleigh Z difference', []};

for m = 1:size(plotmaps,1)
    figure('Position', get(0, 'Screensize'));
    sgtitle([animal ' ' region '-Probe ' ref_condition ' ' plotmaps{m,2} ' (' num2str(ntr_plot) ' trials)'], 'Interpreter', 'none')
    lastax = [];
    for row = 1:nrow
        for q = 1:nshank
            chid = tfmap(row, q);
            ax = subplot(nrow, nshank, (row-1)*nshank + q);
            k = find(chansuse == chid);
            if isempty(k)
                title(['CH' num2str(chid) ' removed'])
                axis off
                continue
            end
            tf = reshape(plotmaps{m,1}(k, fpi, tpi), numel(fpi), numel(tpi)); %freq x time
            if all(isnan(tf(:)))
                title(['CH' num2str(chid) ' no CSD'])
                axis off
                continue
            end
            domask = ~isempty(plotmaps{m,6});
            if domask
                zz = reshape(plotmaps{m,6}(k, fpi, tpi), numel(fpi), numel(tpi));
                tf(zz < rayZ_thresh(1)) = NaN; %blank everything outside the p~0.05 outline
            end
            if ~all(isnan(tf(:))) %contour warns on an all-NaN panel
                contour(times_ms(tpi), frex(fpi), tf, tfplot_levels, 'Fill', 'on')
            end
            set(ax, 'YDir', 'normal', 'CLim', plotmaps{m,4}, ...
                'XLim', times_ms(tpi([1 end])), 'YLim', frex(fpi([1 end])), 'Box', 'on')
            colormap(ax, plotmaps{m,3})
            if domask
                hold on
                for lv = 1:numel(rayZ_thresh) %significance outlines from Rayleigh Z
                    if any(zz(:) >= rayZ_thresh(lv))
                        contour(times_ms(tpi), frex(fpi), zz, rayZ_thresh(lv)*[1 1], 'LineColor', rayZ_colors{lv}, 'LineWidth', 1);
                    end
                end
                hold off
            end
            xline(0, 'LineWidth', 1, 'Color', 'r'); %kept even where the panel is blanked
            title(['CH' num2str(chid)])
            if row < nrow
                set(ax, 'XTickLabel', []) %tick labels on the bottom row only, so they don't hit the titles below
            end
            if q > 1
                set(ax, 'YTickLabel', [])
            end
            if row == nrow && q == 1
                xlabel('Time (ms)')
                ylabel('Frequency (Hz)')
            end
            lastax = ax;
        end
    end
    if ~isempty(lastax)
        axpos = lastax.Position;
        c = colorbar(lastax, 'eastoutside');
        c.Label.String = plotmaps{m,5};
        lastax.Position = axpos; %keep the panel the same size as the others
    end
end

%}

%% 8 (OFF) probe 1 vs probe 2 phase lag, PLV and PSI (volume-conduction check)
%switched off 2026-11-02 with the other test stuff (it also needs the CSD). Kept for later
%{
%Every channel of this probe (probe 1, region) against every channel of a second probe (probe 2,
%p2_region), in time-frequency windows: trial-mean of exp(i*(phase1 - phase2)) gives the PLV
%(consistency of the lag, 0-1) and the mean phase lag. Volume conduction = high PLV at lag ~0 (or ~pi),
%and it should survive ERP removal in the LFP but drop in the CSD. Done for LFP-LFP and CSD-CSD, with
%the ERP kept and removed (plain trial mean subtracted per channel). With the ERP kept, two areas that
%both respond to the stimulus show PLV ~ ITPC1*ITPC2 even with no interaction.
%Also the phase-slope index (PSI; Nolte et al. 2008, Phys Rev Lett 100:234101): the change of coherency
%phase across the band's frequencies. Zero-lag (volume-conducted) mixing has no phase slope, so PSI ~0
%for it; a real delay gives a slope whose sign shows who leads. Standardised by a jackknife over trials.
p2_region = 'S1'; %ProbeInfo.Areas entry of the second probe (source region)
p2_windows = {'alpha', [8 12], [20 100]; ... %name, band (Hz, integer steps within frex), time window (ms re: stim)
    'theta', [4 8], [0 200];
    'lologam', [20 40], [20 40]};
p2_psi_clim = [-4 4]; %colour limits for standardised PSI (|z| > 2 ~ p < 0.05, uncorrected)
p2_lag_zero_deg = 3; %lag panels: cells with |lag| <= this (deg) drawn black (zero lag, volume-conduction-like)
p2_psi_zero = 0.2; %PSI panels: cells with |PSI z| <= this drawn black
%over what stretch of time should I track the coupling moment by moment? (ms re: stim, must contain every window)
p2_tr_ms = [-200 500]; %brief coupling gets diluted in a window average, so this keeps PLV/lag at every time point too
%probe 2 channels (probe channels minus its Ch_Remove), ordered by shank then depth for plotting
prb2 = find(ismember(ProbeInfo.Areas, p2_region));
if numel(prb2) ~= 1 || prb2 == prb
    error('p2_region ''%s'' must match exactly one probe other than %s. ProbeInfo.Areas: %s', ...
        p2_region, region, strjoin(string(ProbeInfo.Areas), ', '));
end
chans2 = ProbeInfo.ProbeIds{prb2}(:);
if numel(ProbeInfo.Ch_Remove) >= prb2
    chans2 = chans2(~ismember(chans2, ProbeInfo.Ch_Remove{prb2}));
end
nch1 = numel(chansuse);
nch2 = numel(chans2);
[row2, shank2] = deal(zeros(nch2,1));
for ch = 1:nch2
    [row2(ch), shank2(ch)] = find(ProbeInfo.ProbeMaps{prb2} == chans2(ch));
end
[row1, shank1] = deal(zeros(nch1,1));
for ch = 1:nch1
    [row1(ch), shank1(ch)] = find(ProbeInfo.ProbeMaps{prb} == chansuse(ch));
end
[~, ord1] = sortrows([shank1 row1]);
[~, ord2] = sortrows([shank2 row2]);

%probe 2 LFP for exactly the trials in trial_table (stim_lfp_stimchunks holds every channel of the recording)
ntr = height(trial_table);
lfp2 = zeros(ntr, nseg, nch2);
cf = unique(trial_table(:, {'condition','file'}), 'rows');
for i = 1:height(cf)
    rows = find(trial_table.condition == cf.condition(i) & trial_table.file == cf.file(i));
    fpath = condi_dir{cf.condition(i)}{cf.file(i)};
    S = load(fpath, 'stim_lfp_stimchunks');
    if size(S.stim_lfp_stimchunks,3) < max(chans2)
        error('%s: stim_lfp_stimchunks has %d channels, probe 2 needs up to %d', fpath, size(S.stim_lfp_stimchunks,3), max(chans2));
    end
    lfp2(rows,:,:) = double(S.stim_lfp_stimchunks(trial_table.trial(rows), seg_idx, chans2));
    clear S
end
csd2 = CSD_flex(lfp2, ProbeInfo, prb2, 'ChRemove', true, ...
    'spacing', csd_spacing_um, 'radius', csd_radius_um, 'conductivity', csd_conductivity, 'minchans', csd_min_chans);

%wavelets as in section 5 (same cycles per frequency: logspace of range_cycles over frex)
cyc_all = logspace(log10(range_cycles(1)), log10(range_cycles(end)), numel(frex));
wavtime = -2:1/1000:2;
half_wave = (length(wavtime)-1)/2;
nConv = length(wavtime) + nseg*ntr - 1;

nwin = size(p2_windows,1);
p12 = struct();
p12.dims = {'probe1 chan','probe2 chan','window','signal {LFP,CSD}','ERP {kept,removed}'};
p12.plv = nan(nch1, nch2, nwin, 2, 2, 'single');
p12.lag = nan(nch1, nch2, nwin, 2, 2, 'single');
p12.lag_ms = nan(nch1, nch2, nwin, 2, 2, 'single');
p12.psi = nan(nch1, nch2, nwin, 2, 2, 'single');
p12.psi_z = nan(nch1, nch2, nwin, 2, 2, 'single');

%time-resolved version: same PLV/lag but kept at every time point in p2_tr_ms (averaged over the band only)
tridx = find(seg_ms >= p2_tr_ms(1) & seg_ms <= p2_tr_ms(2)); %positions of p2_tr_ms inside the padded segment
ntt = numel(tridx);
p12.dims_t = {'probe1 chan','probe2 chan','time','window','signal {LFP,CSD}','ERP {kept,removed}'};
p12.times_tr = seg_ms(tridx); %ms re: stim
p12.plv_t = nan(nch1, nch2, ntt, nwin, 2, 2, 'single'); % [probe1 chans x probe2 chans x time x windows x signals x ERP modes]
p12.lag_t = nan(nch1, nch2, ntt, nwin, 2, 2, 'single');
for w = 1:nwin
    fw = p2_windows{w,2}(1):p2_windows{w,2}(2);
    nfw = numel(fw);
    if nfw < 2
        error('p2_windows ''%s'': PSI needs a band with at least 2 frequencies (1 Hz steps).', p2_windows{w,1});
    end
    if p2_windows{w,3}(1) < p2_tr_ms(1) || p2_windows{w,3}(2) > p2_tr_ms(2)
        error('p2_windows ''%s'' time window must sit inside p2_tr_ms [%d %d].', p2_windows{w,1}, p2_tr_ms);
    end
    iswin = seg_ms(tridx) >= p2_windows{w,3}(1) & seg_ms(tridx) <= p2_windows{w,3}(2); %which tracked time points fall in this window
    ntw = sum(iswin);
    for sg = 1:2 %1 = LFP-LFP, 2 = CSD-CSD
        if sg == 1
            X1 = lfp_ref; X2 = lfp2;
        else
            X1 = csd_ref; X2 = csd2;
        end
        Csum = zeros(nch1, nch2, 2); %sum over freqs/times of the trial-mean vector, per ERP mode
        Psum = zeros(nch1, nch2, 2); %sum of its length (PLV)
        Sk = complex(zeros(nch1, nch2, ntr, nfw, 2)); %per-trial cross-spectrum (summed over window times), for PSI
        P1k = zeros(nch1, ntr, nfw, 2); %per-trial autospectra
        P2k = zeros(nch2, ntr, nfw, 2);
        Ctsum = complex(zeros(nch1, nch2, ntt, 2)); %sum over freqs of the trial-mean vector at each time point
        Ptsum = zeros(nch1, nch2, ntt, 2); %sum over freqs of its length
        for fi = 1:nfw
            f = fw(fi);
            sf = interp1(frex, cyc_all, f) / (2*pi*f);
            wX = fft(exp(2*1i*pi*f.*wavtime) .* exp(-wavtime.^2./(2*sf^2)), nConv);
            wX = wX ./ max(abs(wX));
            %complex wavelet coefficients over p2_tr_ms (the window is a piece of it), trials x time x channels, for both probes
            A1 = nan(ntr, ntt, nch1);
            A2 = nan(ntr, ntt, nch2);
            for pr = 1:2
                if pr == 1, X = X1; nc = nch1; else, X = X2; nc = nch2; end
                for ch = 1:nc
                    dat = X(:,:,ch);
                    if all(isnan(dat(:)))
                        continue %shank with no CSD
                    end
                    dat = dat - mean(dat, 2); %per-trial demean, as in section 5
                    as = ifft(wX .* fft(reshape(dat', 1, []), nConv));
                    as = reshape(as(half_wave+1:end-half_wave), nseg, ntr);
                    if pr == 1, A1(:,:,ch) = as(tridx,:).'; else, A2(:,:,ch) = as(tridx,:).'; end
                end
            end
            for e = 1:2
                if e == 1
                    B1 = A1; B2 = A2;
                else %ERP removed: convolution is linear, so subtracting the trial-mean coefficients = subtracting the ERP before the wavelet
                    B1 = A1 - mean(A1, 1); B2 = A2 - mean(A2, 1);
                end
                U1 = B1 ./ abs(B1); U2 = B2 ./ abs(B2);
                %PSI ingredients: per-trial cross-spectrum and autospectra, summed over the window's time points
                Sk(:,:,:,fi,e) = pagemtimes(permute(B1(:,iswin,:), [3 2 1]), permute(conj(B2(:,iswin,:)), [2 3 1])); %nch1 x nch2 x ntr
                P1k(:,:,fi,e) = reshape(sum(abs(B1(:,iswin,:)).^2, 2), ntr, nch1).';
                P2k(:,:,fi,e) = reshape(sum(abs(B2(:,iswin,:)).^2, 2), ntr, nch2).';
                for t = 1:ntt
                    C = reshape(U1(:,t,:), ntr, nch1).' * conj(reshape(U2(:,t,:), ntr, nch2)) / ntr; %nch1 x nch2, mean of exp(i(phi1-phi2))
                    Ctsum(:,:,t,e) = Ctsum(:,:,t,e) + C; %time-resolved: keep every time point
                    Ptsum(:,:,t,e) = Ptsum(:,:,t,e) + abs(C);
                    if iswin(t) == 1 %window average: only the window's time points
                        Csum(:,:,e) = Csum(:,:,e) + C;
                        Psum(:,:,e) = Psum(:,:,e) + abs(C);
                    end
                end
            end
        end
        for e = 1:2
            p12.plv(:,:,w,sg,e) = single(Psum(:,:,e) / (numel(fw)*ntw));
            p12.lag(:,:,w,sg,e) = single(angle(Csum(:,:,e))); %positive = probe 1 phase ahead of probe 2
            p12.lag_ms(:,:,w,sg,e) = single(angle(Csum(:,:,e)) / (2*pi*mean(fw)) * 1000);
            p12.plv_t(:,:,:,w,sg,e) = single(Ptsum(:,:,:,e) / nfw); %band-averaged PLV at each time point
            p12.lag_t(:,:,:,w,sg,e) = single(angle(Ctsum(:,:,:,e))); %band-averaged lag at each time point, radians

            %PSI = imag(sum over adjacent frequencies of conj(coherency(f)) * coherency(f+1)); coherency
            %phase = phi1 - phi2, so PSI > 0 = probe 1 leads. Jackknife (leave one trial out) for its SD.
            S = sum(Sk(:,:,:,:,e), 3); %nch1 x nch2 x 1 x nfw
            P1 = sum(P1k(:,:,:,e), 2); %nch1 x 1 x nfw
            P2 = sum(P2k(:,:,:,e), 2);
            psi = zeros(nch1, nch2);
            psi_loo = zeros(nch1, nch2, ntr);
            for fi = 1:nfw
                Cf = S(:,:,1,fi) ./ sqrt(P1(:,1,fi) .* P2(:,1,fi).');
                Cl = (S(:,:,1,fi) - Sk(:,:,:,fi,e)) ./ sqrt(reshape(P1(:,1,fi) - P1k(:,:,fi,e), nch1, 1, ntr) ...
                    .* reshape(P2(:,1,fi) - P2k(:,:,fi,e), 1, nch2, ntr));
                if fi > 1
                    psi = psi + imag(conj(Cprev) .* Cf);
                    psi_loo = psi_loo + imag(conj(Clprev) .* Cl);
                end
                Cprev = Cf;
                Clprev = Cl;
            end
            psi_sd = sqrt((ntr-1)/ntr * sum((psi_loo - mean(psi_loo, 3)).^2, 3));
            p12.psi(:,:,w,sg,e) = single(psi);
            p12.psi_z(:,:,w,sg,e) = single(psi ./ psi_sd);
        end
    end
end
p12.chans1 = chansuse;
p12.chans2 = chans2;
p12.region1 = region;
p12.region2 = p2_region;
p12.windows = p2_windows;
p12.ntrials = ntr;
clear X X1 X2 A1 A2 B1 B2 U1 U2 C Csum Psum Ctsum Ptsum dat as wX Sk P1k P2k S P1 P2 Cf Cl Cprev Clprev psi psi_loo psi_sd

%plots: one figure per window; rows LFP / CSD, columns PLV, lag and PSI (z) with the ERP kept / removed.
%Lag is faded where PLV is low (AlphaData = PLV / max PLV of that panel)
sglab = {'LFP-LFP', 'CSD-CSD'};
erplab = {'ERP kept', 'ERP removed'};
b1 = find(diff(shank1(ord1))) + 0.5; %shank boundaries in plot order
b2 = find(diff(shank2(ord2))) + 0.5;
for w = 1:nwin
    fc = mean(p2_windows{w,2}(1):p2_windows{w,2}(2));
    figure('Position', get(0, 'Screensize'));
    sgtitle(sprintf('%s %s vs %s, %s %d-%d Hz, %d to %d ms (%d trials). Lag + = %s phase ahead; faded where PLV is low', ...
        animal, region, p2_region, p2_windows{w,1}, p2_windows{w,2}, p2_windows{w,3}, ntr, region), 'Interpreter', 'none')
    for sg = 1:2
        for e = 1:2
            plv = double(p12.plv(ord1, ord2, w, sg, e));
            lagd = rad2deg(double(p12.lag(ord1, ord2, w, sg, e)));
            psiz = double(p12.psi_z(ord1, ord2, w, sg, e));
            for kind = 1:3 %1 = PLV, 2 = lag, 3 = PSI
                ax = subplot(2, 6, (sg-1)*6 + (e-1)*3 + kind);
                if kind == 1
                    imagesc(plv, 'AlphaData', ~isnan(plv));
                    set(ax, 'CLim', [0 1]);
                    colormap(ax, parula(256));
                    c = colorbar(ax);
                    c.Label.String = 'PLV';
                    title([sglab{sg} ' PLV, ' erplab{e}])
                elseif kind == 3
                    imagesc(psiz, 'AlphaData', ~isnan(psiz));
                    set(ax, 'CLim', p2_psi_clim);
                    colormap(ax, slanCM('coolwarm'));
                    c = colorbar(ax);
                    c.Label.String = sprintf('PSI (z); black = |z| <= %g', p2_psi_zero);
                    title(sprintf('%s PSI, %s (+ = %s leads)', sglab{sg}, erplab{e}, region), 'Interpreter', 'none')
                    hold on %black layer over near-zero PSI
                    image(zeros([size(psiz) 3]), 'AlphaData', double(abs(psiz) <= p2_psi_zero));
                    hold off
                else
                    lagalpha = fillmissing(plv / max(plv(:), [], 'omitnan'), 'constant', 0);
                    imagesc(lagd, 'AlphaData', lagalpha);
                    set(ax, 'CLim', [-180 180]);
                    colormap(ax, hsv(256)); %circular: -180 and 180 deg are the same colour
                    c = colorbar(ax);
                    c.Label.String = sprintf('lag (deg); black = |lag| <= %g deg', p2_lag_zero_deg);
                    title(sprintf('%s lag, %s (90 deg = %.1f ms)', sglab{sg}, erplab{e}, 250/fc))
                    hold on %black layer over near-zero lag, faded by PLV like the lag itself
                    image(zeros([size(lagd) 3]), 'AlphaData', double(abs(lagd) <= p2_lag_zero_deg) .* lagalpha);
                    hold off
                end
                hold on
                for bb = b1', yline(bb, 'k', 'LineWidth', 1); end
                for bb = b2', xline(bb, 'k', 'LineWidth', 1); end
                hold off
                set(ax, 'YTick', 1:nch1, 'YTickLabel', chansuse(ord1), 'XTick', 1:nch2, 'XTickLabel', chans2(ord2), 'FontSize', 6)
                ylabel([region ' channel'])
                xlabel([p2_region ' channel'])
            end
        end
    end
end

%time course of the coupling: one figure per window, rows LFP / CSD, columns ERP kept / removed.
%every channel pair gets its own thin see-through line, so brief coupling in just a few pairs still shows up
plv_chance = 0.886/sqrt(ntr); %roughly what PLV looks like with no coupling at all, for this many trials
for w = 1:nwin
    figure();
    set(gcf, 'Position', get(0, 'Screensize'));
    sgtitle(sprintf('%s %s vs %s, %s %d-%d Hz: PLV over time, all channel pairs (%d trials). Grey = analysis window, dashed = chance', ...
        animal, region, p2_region, p2_windows{w,1}, p2_windows{w,2}, ntr), 'Interpreter', 'none')
    for sg = 1:2
        for e = 1:2
            subplot(2, 2, (sg-1)*2 + e)
            pt = reshape(double(p12.plv_t(:,:,:,w,sg,e)), nch1*nch2, ntt); % [channel pairs x time]
            pt = pt(~all(isnan(pt), 2), :); %drop pairs with no data (no-CSD shanks)
            %old version: one median line + middle-50% band. Replaced with every pair, because a median hides coupling in only a few pairs
            %pmed = median(pt, 1, 'omitnan');
            %plo = prctile(pt, 25, 1); %prctile skips NaN pairs (no-CSD shanks)
            %phi = prctile(pt, 75, 1);
            linealpha = max(0.03, min(0.5, 10/size(pt,1))); %more pairs = fainter lines, so dense bundles show where most pairs sit
            hold on
            patch([p2_windows{w,3}(1) p2_windows{w,3}(2) p2_windows{w,3}(2) p2_windows{w,3}(1)], [0 0 1 1], [.9 .9 .9], 'EdgeColor', 'none'); %analysis window
            %fill([p12.times_tr fliplr(p12.times_tr)], [plo fliplr(phi)], [.6 .75 1], 'EdgeColor', 'none'); %middle 50% of pairs
            %plot(p12.times_tr, pmed, 'b', 'LineWidth', 1.5); %median pair
            plot(p12.times_tr, pt.', 'Color', [0 0 1 linealpha], 'LineWidth', 0.75); %one line per channel pair
            yline(plv_chance, '--k');
            xline(0, 'LineWidth', 1, 'Color', 'r');
            hold off
            xlim(p2_tr_ms)
            ylim([0 1])
            box on
            xlabel('Time (ms)')
            ylabel('PLV')
            title(sprintf('%s, %s (%d pairs)', sglab{sg}, erplab{e}, size(pt,1)))
        end
    end
end
%}



%Now beginning the plotting section

%% 9 Data Visualization: pre vs post power spectra in the probe layout
%black = pre-stim, red = post-stim. One figure per probe per window size
disp('Plotting: power spectra')
speclist = {spec_pre500, spec_post500, 500; ...
    spec_pre200, spec_post200, 200;
    spec_pre100, spec_post100, 100}; %pre, post, window size (ms)
for prb = prbs
    probemap = ProbeInfo.ProbeMaps{prb}; % [depth rows x shanks]
    [nrow, nshank] = size(probemap);
    onprobe = ch_probe == prb; %this probe's channels
    for w = 1:size(speclist,1)
        spre = speclist{w,1};
        spost = speclist{w,2};
        both = [spre(onprobe,:) spost(onprobe,:)];
        yl = [min(both, [], 'all') max(both, [], 'all')]; %same y range for every channel on the probe

        figure();
        set(gcf, 'Position', get(0, 'Screensize'));
        sgtitle([animal ' ' ProbeInfo.Areas{prb} '-Probe ' ref_condition ' power spectrum: -' num2str(speclist{w,3}) ' to 0 ms (black) vs 0 to ' ...
            num2str(speclist{w,3}) ' ms (red), ' num2str(ntr) ' trials'], 'Interpreter', 'none')
        firstpanel = 1;
        for row = 1:nrow
            for shank = 1:nshank
                chid = probemap(row, shank);
                subplot(nrow, nshank, (row-1)*nshank + shank)
                k = find(chansuse == chid & ch_probe == prb); %where this channel lives in the channel dim
                if isempty(k) == true
                    title(['CH' num2str(chid) ' removed'])
                    axis off
                    continue
                end
                semilogy(frex, spre(k,:), 'k', 'LineWidth', 1) %pre
                hold on
                semilogy(frex, spost(k,:), 'r', 'LineWidth', 1) %post
                hold off
                xlim([frex(1) frex(end)])
                ylim(yl)
                title(['CH' num2str(chid)])
                if row < nrow
                    set(gca, 'XTickLabel', []) %tick labels on the bottom row only, so they don't hit the titles below
                end
                if row == nrow && shank == 1
                    xlabel('Frequency (Hz)')
                    ylabel('Power (uV^2)')
                end
                if firstpanel == 1
                    legend({'pre', 'post'}, 'Location', 'northeast')
                    firstpanel = 0;
                end
            end
        end
    end
end

%% 10 Data Visualization: ITPC in the probe layout, with Rayleigh outlines
%cyan outline = Rayleigh p ~ 0.05, white = p ~ 0.01 (per point, uncorrected). Thresholds come from the trial count
disp('Plotting: ITPC')
tpi = find(times_ms >= itpcplot_ms(1) & times_ms <= itpcplot_ms(2)); %time points to show
fpi = find(frex >= itpcplot_hz(1) & frex <= itpcplot_hz(2)); %frequencies to show

%Rayleigh Z at p = 0.05 and 0.01 for this trial count: the Zar formula from section 6,
%p = exp(sqrt(1 + 4n + 4(n^2 - Rn^2)) - (1 + 2n)), solved for Z = Rn^2/n
rayZ_p = [0.05 0.01];
L = log(rayZ_p) + 1 + 2*ntr;
rayZ_thresh = (1 + 4*ntr + 4*ntr^2 - L.^2) / (4*ntr);
rayZ_colors = {'c', 'w'}; %outline colour for each p level

for prb = prbs
    probemap = ProbeInfo.ProbeMaps{prb}; % [depth rows x shanks]
    [nrow, nshank] = size(probemap);

    figure();
    set(gcf, 'Position', get(0, 'Screensize'));
    sgtitle(sprintf('%s %s-Probe %s ITPC (%d trials). Cyan p~0.05 (Z=%.2f), white p~0.01 (Z=%.2f), uncorrected', ...
        animal, ProbeInfo.Areas{prb}, ref_condition, ntr, rayZ_thresh), 'Interpreter', 'none')
    lastax = [];
    for row = 1:nrow
        for shank = 1:nshank
            chid = probemap(row, shank);
            ax = subplot(nrow, nshank, (row-1)*nshank + shank);
            k = find(chansuse == chid & ch_probe == prb);
            if isempty(k) == true
                title(['CH' num2str(chid) ' removed'])
                axis off
                continue
            end
            tf = reshape(itpc_lfp(k, fpi, tpi), numel(fpi), numel(tpi)); % [frequencies x time]
            zz = reshape(rayleighZ_lfp(k, fpi, tpi), numel(fpi), numel(tpi)); %Rayleigh Z for the outlines
            contour(times_ms(tpi), frex(fpi), tf, 50, 'Fill', 'on')
            set(ax, 'YDir', 'normal', 'CLim', itpc_clim)
            colormap(ax, parula(256))
            hold on
            for lv = 1:numel(rayZ_thresh) %significance outlines
                if any(zz(:) >= rayZ_thresh(lv))
                    contour(times_ms(tpi), frex(fpi), zz, rayZ_thresh(lv)*[1 1], 'LineColor', rayZ_colors{lv}, 'LineWidth', 1);
                end
            end
            xline(0, 'LineWidth', 1, 'Color', 'r'); %stim onset
            hold off
            title(['CH' num2str(chid)])
            if row < nrow
                set(ax, 'XTickLabel', [])
            end
            if shank > 1
                set(ax, 'YTickLabel', [])
            end
            if row == nrow && shank == 1
                xlabel('Time (ms)')
                ylabel('Frequency (Hz)')
            end
            lastax = ax;
        end
    end
    if isempty(lastax) == false
        axpos = lastax.Position;
        c = colorbar(lastax, 'eastoutside');
        c.Label.String = 'ITPC';
        lastax.Position = axpos; %keep the panel the same size as the others
    end
end

%% 10b Data Visualization: significant ITPC only (Rayleigh p < 0.05), the rest blanked
%same layout as section 10, but anything outside the p~0.05 outline is blank, so only the significant
%time-frequency regions are left. Outlines: cyan = p ~ 0.05, white = p ~ 0.01 (per point, uncorrected)
disp('Plotting: significant ITPC only')
tpi = find(times_ms >= itpcplot_ms(1) & times_ms <= itpcplot_ms(2)); %time points to show (rebuilt here so this section runs on its own)
fpi = find(frex >= itpcplot_hz(1) & frex <= itpcplot_hz(2)); %frequencies to show

%Rayleigh Z at p = 0.05 and 0.01 for this trial count (same Zar formula as sections 6 and 10)
rayZ_p = [0.05 0.01];
L = log(rayZ_p) + 1 + 2*ntr;
rayZ_thresh = (1 + 4*ntr + 4*ntr^2 - L.^2) / (4*ntr);
rayZ_colors = {'c', 'w'}; %outline colour for each p level

for prb = prbs
    probemap = ProbeInfo.ProbeMaps{prb}; % [depth rows x shanks]
    [nrow, nshank] = size(probemap);

    figure();
    set(gcf, 'Position', get(0, 'Screensize'));
    sgtitle(sprintf('%s %s-Probe %s ITPC where Rayleigh p < 0.05 (%d trials). Cyan p~0.05 (Z=%.2f), white p~0.01 (Z=%.2f), uncorrected', ...
        animal, ProbeInfo.Areas{prb}, ref_condition, ntr, rayZ_thresh), 'Interpreter', 'none')
    lastax = [];
    for row = 1:nrow
        for shank = 1:nshank
            chid = probemap(row, shank);
            ax = subplot(nrow, nshank, (row-1)*nshank + shank);
            k = find(chansuse == chid & ch_probe == prb);
            if isempty(k) == true
                title(['CH' num2str(chid) ' removed'])
                axis off
                continue
            end
            tf = reshape(itpc_lfp(k, fpi, tpi), numel(fpi), numel(tpi)); % [frequencies x time]
            zz = reshape(rayleighZ_lfp(k, fpi, tpi), numel(fpi), numel(tpi)); %Rayleigh Z, same shape
            tf(zz < rayZ_thresh(1)) = NaN; %blank everything outside the p~0.05 outline

            if any(isnan(tf(:)) == false) %only fill if something survived (contour complains about all-NaN panels)
                contour(times_ms(tpi), frex(fpi), tf, 50, 'Fill', 'on')
            end
            set(ax, 'YDir', 'normal', 'CLim', itpc_clim, ... %same colour scale as section 10
                'XLim', times_ms(tpi([1 end])), 'YLim', frex(fpi([1 end])), 'Box', 'on') %explicit limits, so blank panels keep their axes
            colormap(ax, parula(256))
            hold on
            for lv = 1:numel(rayZ_thresh) %significance outlines
                if any(zz(:) >= rayZ_thresh(lv))
                    contour(times_ms(tpi), frex(fpi), zz, rayZ_thresh(lv)*[1 1], 'LineColor', rayZ_colors{lv}, 'LineWidth', 1);
                end
            end
            xline(0, 'LineWidth', 1, 'Color', 'r'); %stim onset, kept even on blank panels
            hold off
            title(['CH' num2str(chid)])
            if row < nrow
                set(ax, 'XTickLabel', [])
            end
            if shank > 1
                set(ax, 'YTickLabel', [])
            end
            if row == nrow && shank == 1
                xlabel('Time (ms)')
                ylabel('Frequency (Hz)')
            end
            lastax = ax;
        end
    end
    if isempty(lastax) == false
        axpos = lastax.Position;
        c = colorbar(lastax, 'eastoutside');
        c.Label.String = 'ITPC';
        lastax.Position = axpos; %keep the panel the same size as the others
    end
end

%% 11 Data Visualization: band-passed ERP on top of the ITPC
%pick a band (probably from the ITPC plots) and its ERP gets drawn over the ITPC, with its own y axis on the right.
%Uses tpi, fpi, rayZ_thresh and rayZ_colors from section 10, so run that first
ishsokay = 0;
while ishsokay == 0
    dlgtitle = 'Band ERP';
    promt = {'Which frequency band? ([low high] in Hz)'};
    fieldsize = [1 50];
    definput = {num2str(erp_band)};
    opts.Resize = 'on';
    opts.WindowStyle = 'normal';
    bandpick = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
    erp_band = str2num(bandpick{1}); %Hz

    %ERP first, then bandpass: the bandpass is a linear zero-phase filter, so filtering the trial average
    %gives the same thing as averaging the filtered trials (and is way faster)
    erp_raw = reshape(mean(lfp_ref, 1), nseg, nch); % [segment x channels], uV
    erp_band_lfp = bandpass(erp_raw, erp_band, 1000, 'ImpulseResponse', 'iir')'; % [channels x segment], uV
    eidx = find(seg_ms >= itpcplot_ms(1) & seg_ms <= itpcplot_ms(2)); %ERP points that line up with the ITPC time axis

    for prb = prbs
        probemap = ProbeInfo.ProbeMaps{prb}; % [depth rows x shanks]
        [nrow, nshank] = size(probemap);
        onprobe = ch_probe == prb;
        elim = max(abs(erp_band_lfp(onprobe, eidx)), [], 'all'); %same ERP scale for the whole probe

        figure();
        set(gcf, 'Position', get(0, 'Screensize'));
        sgtitle(sprintf('%s %s-Probe %s ITPC with %g-%g Hz ERP (black, right axis, uV). Dashed = band edges', ...
            animal, ProbeInfo.Areas{prb}, ref_condition, erp_band), 'Interpreter', 'none')
        for row = 1:nrow
            for shank = 1:nshank
                chid = probemap(row, shank);
                ax = subplot(nrow, nshank, (row-1)*nshank + shank);
                k = find(chansuse == chid & ch_probe == prb);
                if isempty(k) == true
                    title(['CH' num2str(chid) ' removed'])
                    axis off
                    continue
                end
                tf = reshape(itpc_lfp(k, fpi, tpi), numel(fpi), numel(tpi)); % [frequencies x time]
                zz = reshape(rayleighZ_lfp(k, fpi, tpi), numel(fpi), numel(tpi));

                yyaxis left %ITPC on the left axis
                contour(times_ms(tpi), frex(fpi), tf, 50, 'Fill', 'on')
                set(ax, 'CLim', itpc_clim)
                colormap(ax, parula(256))
                hold on
                for lv = 1:numel(rayZ_thresh)
                    if any(zz(:) >= rayZ_thresh(lv))
                        contour(times_ms(tpi), frex(fpi), zz, rayZ_thresh(lv)*[1 1], 'LineColor', rayZ_colors{lv}, 'LineWidth', 1);
                    end
                end
                yline(erp_band(1), '--w'); %band edges
                yline(erp_band(2), '--w');
                hold off
                ylim([frex(fpi(1)) frex(fpi(end))])

                yyaxis right %ERP on the right axis
                plot(seg_ms(eidx), erp_band_lfp(k, eidx), 'k', 'LineWidth', 1)
                if elim > 0
                    ylim([-elim elim])
                end
                ax.YAxis(1).Color = 'k'; %both axes black instead of yyaxis' default blue/orange
                ax.YAxis(2).Color = 'k';
                xline(0, 'LineWidth', 1, 'Color', 'r'); %stim onset
                xlim(itpcplot_ms)

                title(['CH' num2str(chid)])
                if row < nrow
                    set(ax, 'XTickLabel', [])
                end
                if row == nrow && shank == 1
                    xlabel('Time (ms)')
                    yyaxis left
                    ylabel('Frequency (Hz)')
                    yyaxis right
                    ylabel('ERP (uV)')
                end
            end
        end
    end

    %are you happy with that? if yes, continue. If no, pick another band
    quest = 'Happy with that band?';
    btn1 = 'Yes, happy as a clam';
    btn2 = 'No, lemme pick again';
    answer = questdlg(quest, 'Band ERP', btn1, btn2, btn1);
    switch answer
        case btn1
            ishsokay = 1;
        case btn2
            ishsokay = 0;
    end
end

%%
 save('20261011-p7_WactiveTTX_itpc_res.mat', 'erp_raw', 'itpc_lfp', 'itpc2_lfp', 'rayleighZ_lfp')