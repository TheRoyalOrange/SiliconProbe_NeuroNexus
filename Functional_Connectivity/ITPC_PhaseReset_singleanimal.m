% ITPC_PhaseReset_singleanimal.m
%
% Description: For one animal and one probe/region, runs a short-cycle
%   complex Morlet wavelet on peri-stimulus LFP trials from a set of
%   conditions, keeps the per-trial phase of every frequency/time point on
%   every non-Ch_Remove channel, and from it computes:
%     - ITPC and preferred phase angle (prefAngle) per channel x freq x
%       time, per condition group (sec 3), plotted in probe layout (sec 4)
%     - optional ITPC significance (Rayleigh test + BH FDR per channel,
%       sec 5) to help choose time-frequency windows
%     - user-chosen time-frequency windows (sec 6, set AFTER looking at
%       the ITPC plots), per-trial phases inside them (sec 7)
%     - window ITPC as one value per channel, probe/shank plots, optional
%       per-channel permutation test between two groups (sec 8)
%     - a per-trial phase-reset score for EVERY trial (all groups),
%       measured against the reference group's preferred angle in each
%       window (leave-one-out for reference-group trials) (sec 9). The
%       scoring metric is PROVISIONAL (reset_method) - the phase
%       difference itself (reset_d) is stored so the metric can change.
%     - intersite phase clustering (ISPC) between channel pairs inside each
%       window, per group (sec 10)
%   Sections 1-3 must run in order. Then look at sec 4/5 plots, set the
%   windows in sec 6, and run 6 onward.
%   Wavelet method copied from OpenEphys_BaseAnalysis_Bundled.m (TF
%   section: logspace cycles, wavtime -2:2 s, FFT convolution of trials
%   concatenated per channel, waveletX./max(waveletX)), with two
%   differences: range_cycles defaults to [3 5] (BaseAnalysis uses [6 10])
%   for better time resolution, and each trial segment is demeaned before
%   convolution (short-cycle wavelets leak a little DC).
%
% Inputs:
%   ProbeInfo (struct) - loaded from <probeinfo_dir>\<animal>-ProbeInfo.mat
%     .Areas (cell of char) - region per probe; region must match one
%     .ProbeIds (cell, 1 x probenum) - raw channel IDs per probe (these are
%       raw data rows, NOT 1:N per probe - see oshkosh_detect_funct.m)
%     .Ch_Remove (cell, 1 x probenum) - raw channel IDs to exclude
%     .ProbeMaps (cell, 1 x probenum) - 2D [depth rows x shanks] matrix of
%       raw channel IDs. INFERRED: row 1 is one end of the shank (which
%       end - superficial or deep - is not stated anywhere)
%   <animal>-<stim>-LFP.mat files (chosen per condition with uipickfiles),
%     from OpenEphys_BaseAnalysis_Bundled.m:
%     stim_lfp_stimchunks (double, trials x 15000 x channels) - LFP in
%       microvolts at 1 kHz; epoch = stim-5000 : stim+9999, so stimulus
%       onset is sample 5001 (onset_idx). Channel index = raw channel ID.
%     stim_times (double, trials x 1) - onsets in 1 kHz samples
%     tr_keep (int vector) - kept trial indices
%     tr_remove - INFERRED, two formats exist: BaseAnalysis saves removed
%       trial indices (already excluded from tr_keep); ConditionalTrial-
%       Remove_callable.m migrates missing/empty ones to a 0/1 mask the
%       same length as tr_keep. Detected per file (report.files.tr_remove_format).
%     tr_remove_conditional (table, optional) - Name/TrialIdx rows from
%       ConditionalTrialRemove_callable.m; TrialIdx = original trial
%       numbers, NaN = condition checked, nothing removed. Applied only for
%       names listed in tr_conditional_names.
%   User parameters - see section 1.
%
% Outputs (workspace; section 12 saves them):
%   itpc_info (struct) - fields:
%     .ITPC (cell, 1 x ngroups) - each single, nch x nfrex x ntime
%     .prefAngles (cell, 1 x ngroups) - same size, radians
%     .dimfeatures {'chans','freqs','time'}
%     .stimtime (int) - index of t = 0 ms along the time dim
%     .frex (double, 1 x nfrex) - Hz; .times_ms (1 x ntime) - ms re: stim
%     .range_cycles, .windows (struct array: name, freq_hz, time_ms, fidx, tidx)
%     .ch_index (raw channel IDs, nch x 1) - row k of every channel dim
%     .ch_shank, .ch_row (nch x 1) - column/row of each channel in ProbeMaps
%     .itpc_win (cell, 1 x nwin) - nch x ngroups, mean window ITPC
%     .rayleighZ (cell, 1 x nwin) - nch x ngroups, n * itpc_win.^2
%     .itpc_p / .itpc_sig (cell, 1 x ngroups) - nch x nfrex x ntime Rayleigh
%       p and BH-FDR significance ([] if do_signif = 0)
%     .perm_p / .perm_q (cell, 1 x nwin) - nch x 1, window ITPC difference
%       between compare_groups ([] if not run)
%     .trial_table (table) - one row per trial per group: group, condition,
%       file, original trial number, stim time; row order = trial dim of
%       win_u / reset_* arrays
%     .win_u (cell, 1 x nwin) - complex single, nch x nfw x ntw x ntrials,
%       unit phase vectors inside each window
%     .reset_d (cell, 1 x nwin) - single, same size, phase difference
%       trial - reference angle, radians (-pi, pi]
%     .reset_score, .reset_offset (cell, 1 x nwin) - nch x ntrials
%     .reset_method (char)
%     .ispc (cell, 1 x nwin, each cell 1 x ngroups) - nch x nch
%     plus naming fields kept from Plotting_ITPC_probe_singleanimal.m:
%     .condition_directory, .all_conditions, .condition_grouping,
%     .condition_grouping_index, .groupnames, .region, .report
%   phases (cell, 1 x ngroups) - single, nch x nfrex x ntime x ntrials(grp),
%     radians. Workspace only (large), NOT saved.
%   CSVs for R (section 12): <base>_windowITPC.csv (one row per window x
%     group x channel) and <base>_trialReset.csv (one row per window x
%     trial x channel).
%
% Notes:
%   - ITPC is biased upward for small n. Compare groups with rayleighZ or
%     the permutation test (sec 8), not raw ITPC alone.
%   - Neighbouring LFP channels share volume-conducted signal, so ISPC is
%     high near the diagonal regardless of coupling; the between-group
%     difference map is the more meaningful readout.
%   - If one condition is assigned to two groups its trials appear in both;
%     the permutation test and leave-one-out reference then double-count.
%   - Probe layout plots use subplot(depth rows, shanks). The older scripts
%     use the transposed map for the grid, which only matches for square
%     (8x8) maps.
%
% Dependencies: OpenEphys_BaseAnalysis_Bundled.m (LFP files, ProbeInfo),
%   OpenEphys_editProbeInfo_ChRemove.m (Ch_Remove), uipickfiles.m;
%   optional ConditionalTrialRemove_callable.m (tr_remove_conditional).

%% 1 PARAMETERS
animal = '20260820-p8';
region = 'V1'; %must match one entry of ProbeInfo.Areas

condis = {'W', 'W_TTX'}; %condition labels, in the order files will be picked
groupassign = {[1], [2]}; %which conditions make each group
groupnames = {'Whisker', 'Whisker_TTX'};
ref_group = 1; %group whose preferred angle defines "reset"

frex = 2:1:80; %Hz, integer steps so frequency index is unambiguous
range_cycles = [3 5]; %wavelet cycles at lowest/highest frequency
time_ms = [-1000 1000]; %analysis window kept after convolution, ms re: stim
pad_ms = 1000; %extra data convolved on each side, keeps edge effects out of time_ms
onset_idx = 5001; %stimulus onset sample in stim_lfp_stimchunks

tr_conditional_names = {}; %e.g. {'whisker_twitch_artifact'}; extra trial exclusions

do_signif = 1; %run significance tests (sec 5 and sec 8)
alpha = 0.05;
compare_groups = [1 2]; %groups compared in sec 8 (first minus second)
nperm = 1000;

max_phase_gb = 60; %stop before loading if stored phases would exceed this
null_check = 0; %1 = circularly shift each trial by a random lag before the wavelet (sanity check: ITPC should vanish)

probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';
lfp_dir = 'E:\Roy\Processed Silicon Probe Data\LFP';
save_root = 'E:\Roy\Silcon Probe Data Followup Analyses\ITPC';

%% 2 SETUP: channels, files, trials
report = struct();
report.notes = {};

load(fullfile(probeinfo_dir, [animal '-ProbeInfo.mat']), 'ProbeInfo');
prb = find(ismember(ProbeInfo.Areas, region));
if numel(prb) ~= 1
    error('Region %s matches %d probes in ProbeInfo.Areas', region, numel(prb));
end

%channels: raw IDs on this probe minus Ch_Remove. Row k of every channel dim = chansuse(k)
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

%pick LFP files per condition
condi_dir = cell(1,numel(condis));
for con = 1:numel(condis)
    picked = uipickfiles('FilterSpec', fullfile(lfp_dir, animal), 'Prompt', ['Choose ' condis{con} ' LFP files']);
    if ~iscell(picked) || isempty(picked)
        error('No files chosen for condition %s', condis{con});
    end
    condi_dir{con} = picked;
end

%trials used per file: tr_keep, minus tr_remove (either format), minus chosen conditional removals
file_trials = cell(1,numel(condis)); %file_trials{con}{file} = original trial numbers used
file_stimtimes = cell(1,numel(condis));
filerows = {};
for con = 1:numel(condis)
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
                T = C.tr_remove_conditional;
                missing = setdiff(string(tr_conditional_names), string(T.Name));
                if ~isempty(missing)
                    report.notes{end+1} = sprintf('%s: conditional name(s) not recorded: %s', fname, strjoin(missing, ', '));
                end
                rm = T.TrialIdx(ismember(string(T.Name), string(tr_conditional_names)));
                rm = rm(~isnan(rm)); %NaN = checked, nothing removed
                ncond_removed = sum(ismember(tr_use, rm));
                tr_use = tr_use(~ismember(tr_use, rm));
            else
                report.notes{end+1} = sprintf('%s: no tr_remove_conditional table, conditional removals skipped', fname);
            end
        end

        file_trials{con}{file} = tr_use;
        file_stimtimes{con}{file} = F.stim_times(tr_use);
        filerows(end+1,:) = {condis{con}, file, string(fname), numel(F.stim_times), numel(tr_keep), rmformat, ncond_removed, numel(tr_use)}; %#ok<SAGROW>
    end
end
report.files = cell2table(filerows, 'VariableNames', ...
    {'condition','file','filename','ntrials_total','ntrials_keep','tr_remove_format','ncond_removed','ntrials_used'});

%trial bookkeeping: one row per trial per group. Row order = trial dim of win_u / reset_* arrays
ngrp = numel(groupassign);
ntr = zeros(1,ngrp);
tt = {};
for grp = 1:ngrp
    rig = 0;
    for c = groupassign{grp}
        for file = 1:numel(condi_dir{c})
            tr = file_trials{c}{file}(:);
            n = numel(tr);
            [~, fname] = fileparts(condi_dir{c}{file});
            tt{end+1} = table(repmat(grp,n,1), repmat(string(groupnames{grp}),n,1), repmat(c,n,1), ...
                repmat(string(condis{c}),n,1), repmat(file,n,1), repmat(string(fname),n,1), tr, ...
                file_stimtimes{c}{file}(:), rig+(1:n)', ...
                'VariableNames', {'group','groupname','condition','condiname','file','filename','trial','stim_time','row_in_group'}); %#ok<SAGROW>
            rig = rig + n;
        end
    end
    ntr(grp) = rig;
end
trial_table = vertcat(tt{:});
report.ntrials_per_group = ntr;
if numel([groupassign{:}]) ~= numel(unique([groupassign{:}]))
    report.notes{end+1} = 'A condition is assigned to more than one group: its trials are double-counted in the permutation test and leave-one-out reference.';
end

%% 3 WAVELET + PHASE (convolution per channel; phases kept for all channels)
nfrex = numel(frex);
times_ms = time_ms(1):time_ms(2);
ntime = numel(times_ms);
seg_idx = onset_idx + ((time_ms(1)-pad_ms):(time_ms(2)+pad_ms)); %samples loaded per trial
nseg = numel(seg_idx);
keep_idx = pad_ms + (1:ntime); %positions of time_ms inside the segment

%wavelet parameters (as in OpenEphys_BaseAnalysis_Bundled.m TF section)
s = logspace(log10(range_cycles(1)), log10(range_cycles(end)), nfrex) ./ (2*pi*frex);
wavtime = -2:1/1000:2;
half_wave = (length(wavtime)-1)/2;
nWave = length(wavtime);

report.phase_gb = 4*nch*nfrex*ntime*sum(ntr)/1e9;
if report.phase_gb > max_phase_gb
    error('Stored phases would need %.1f GB (> max_phase_gb = %g). Narrow frex or time_ms.', report.phase_gb, max_phase_gb);
end
if seg_idx(1) < 1
    error('time_ms(1) - pad_ms reaches before the start of the epoch.');
end
if null_check == 1
    rng(0);
    report.notes{end+1} = 'NULL CHECK ON: each trial was circularly shifted by a random lag before the wavelet. Results are a null, not data.';
end

phases = cell(1,ngrp);
for grp = 1:ngrp
    phases{grp} = zeros(nch, nfrex, ntime, ntr(grp), 'single');
end

for grp = 1:ngrp
    for c = groupassign{grp}
        for file = 1:numel(condi_dir{c})
            disp(['Wavelet: group ' groupnames{grp} ', ' condis{c} ' file ' num2str(file)])
            rows = trial_table.row_in_group(trial_table.group == grp & trial_table.condition == c & trial_table.file == file);
            tr = file_trials{c}{file};
            ntrf = numel(tr);

            S = load(condi_dir{c}{file}, 'stim_lfp_stimchunks'); %not -v7.3, so load once and loop channels
            if size(S.stim_lfp_stimchunks,2) < seg_idx(end) || size(S.stim_lfp_stimchunks,3) < max(chansuse)
                error('%s: stim_lfp_stimchunks is %s, too small for this segment/channel set', condi_dir{c}{file}, mat2str(size(S.stim_lfp_stimchunks)));
            end
            lfp = double(S.stim_lfp_stimchunks(tr, seg_idx, chansuse)); %trials x nseg x nch
            clear S

            %wavelet FFTs for this file's convolution length
            nConv = nWave + nseg*ntrf - 1;
            waveletX = zeros(nfrex, nConv);
            for fi = 1:nfrex
                wavelet = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
                wX = fft(wavelet, nConv);
                waveletX(fi,:) = wX ./ max(wX);
            end

            for ch = 1:nch
                dat = lfp(:,:,ch); %trials x nseg
                dat = dat - mean(dat,2); %demean each trial segment
                if null_check == 1
                    for t = 1:ntrf
                        dat(t,:) = circshift(dat(t,:), randi(nseg));
                    end
                end
                dataX = fft(reshape(dat',1,[]), nConv); %trials concatenated end to end

                chph = zeros(nfrex, ntime, ntrf, 'single');
                for fi = 1:nfrex
                    as = ifft(waveletX(fi,:) .* dataX);
                    as = as(half_wave+1:end-half_wave);
                    as = reshape(as, nseg, ntrf); %time x trials
                    chph(fi,:,:) = single(angle(as(keep_idx,:)));
                end
                phases{grp}(ch,:,:,rows) = reshape(chph, [1 nfrex ntime ntrf]);
            end
        end
    end
end
clear lfp waveletX dataX chph as

%ITPC and preferred angle per group
itpc = cell(1,ngrp);
prefAngle = cell(1,ngrp);
for grp = 1:ngrp
    itpc{grp} = zeros(nch, nfrex, ntime, 'single');
    prefAngle{grp} = zeros(nch, nfrex, ntime, 'single');
    for ch = 1:nch
        m = mean(exp(1i*reshape(phases{grp}(ch,:,:,:), nfrex, ntime, [])), 3);
        itpc{grp}(ch,:,:) = abs(m);
        prefAngle{grp}(ch,:,:) = angle(m);
    end
end

%% 4 PLOT: probe-wide TF ITPC per group (look here to choose windows)
plot_time_ms = [-100 300];
plot_freq_hz = [frex(1) frex(end)];
clims = [0 0.5];

probemap = ProbeInfo.ProbeMaps{prb};
[nrow, nshank] = size(probemap);
for grp = 1:numel(itpc)
    figure('Position', get(0,'Screensize'));
    sgtitle([animal ' ' region ' ' groupnames{grp} ' ITPC (n = ' num2str(ntr(grp)) ')'], 'Interpreter', 'none')
    for r = 1:nrow
        for sh = 1:nshank
            ch = find(chansuse == probemap(r,sh));
            if isempty(ch); continue; end
            subplot(nrow, nshank, (r-1)*nshank + sh)
            contourf(times_ms, frex, double(squeeze(itpc{grp}(ch,:,:))), 25, 'LineColor', 'none')
            set(gca, 'YDir', 'normal', 'XLim', plot_time_ms, 'YLim', plot_freq_hz, 'CLim', clims)
            xline(0, 'r', 'LineWidth', 1)
            title(['CH' num2str(probemap(r,sh))])
        end
    end
end

%group difference (compare_groups(1) - compare_groups(2)). Raw ITPC: biased by unequal n
if numel(itpc) >= 2
    g1 = compare_groups(1); g2 = compare_groups(2);
    figure('Position', get(0,'Screensize'));
    sgtitle([animal ' ' region ' ITPC ' groupnames{g1} ' minus ' groupnames{g2}], 'Interpreter', 'none')
    for r = 1:nrow
        for sh = 1:nshank
            ch = find(chansuse == probemap(r,sh));
            if isempty(ch); continue; end
            subplot(nrow, nshank, (r-1)*nshank + sh)
            contourf(times_ms, frex, double(squeeze(itpc{g1}(ch,:,:) - itpc{g2}(ch,:,:))), 25, 'LineColor', 'none')
            set(gca, 'YDir', 'normal', 'XLim', plot_time_ms, 'YLim', plot_freq_hz, 'CLim', [-0.3 0.3])
            xline(0, 'k', 'LineWidth', 1)
            title(['CH' num2str(probemap(r,sh))])
        end
    end
end

%% 4b PLOT: single-channel close-up, all groups
chcheck = chansuse(1); %raw channel ID
plot_time_ms = [-100 300];
plot_freq_hz = [frex(1) frex(end)];
clims = [0 0.5];

ch = find(chansuse == chcheck);
figure();
for grp = 1:numel(itpc)
    subplot(1, numel(itpc), grp)
    contourf(times_ms, frex, double(squeeze(itpc{grp}(ch,:,:))), 40, 'LineColor', 'none')
    set(gca, 'YDir', 'normal', 'XLim', plot_time_ms, 'YLim', plot_freq_hz, 'CLim', clims)
    xline(0, 'r', 'LineWidth', 1)
    colorbar
    title([groupnames{grp} ' CH' num2str(chcheck)], 'Interpreter', 'none')
    xlabel('Time (ms)'); ylabel('Frequency (Hz)')
end

%% 5 ITPC SIGNIFICANCE (optional): Rayleigh p per pixel, BH FDR per channel
itpc_p = cell(1,numel(itpc));
itpc_sig = cell(1,numel(itpc));
if do_signif == 1
    report.sig_channels = zeros(1,numel(itpc)); %channels with any significant post-stim pixel
    for grp = 1:numel(itpc)
        n = ntr(grp);
        R = double(itpc{grp}) * n;
        p = exp(sqrt(1 + 4*n + 4*(n^2 - R.^2)) - (1 + 2*n)); %Rayleigh p (Zar approximation)
        p = min(max(p, 0), 1);
        sig = false(size(p));
        for ch = 1:size(p,1)
            pc = p(ch,:,:);
            [ps, si] = sort(pc(:));
            m = numel(ps);
            q = ps .* m ./ (1:m)';
            q = flipud(cummin(flipud(q)));
            qq = zeros(m,1);
            qq(si) = min(q, 1);
            sig(ch,:,:) = reshape(qq <= alpha, [1 size(p,2) size(p,3)]);
        end
        itpc_p{grp} = single(p);
        itpc_sig{grp} = sig;
        report.sig_channels(grp) = sum(any(reshape(sig(:,:,times_ms > 0), size(sig,1), []), 2));
    end

    %probe layout with significant regions outlined
    plot_time_ms = [-100 300];
    plot_freq_hz = [frex(1) frex(end)];
    clims = [0 0.5];
    probemap = ProbeInfo.ProbeMaps{prb};
    [nrow, nshank] = size(probemap);
    for grp = 1:numel(itpc)
        figure('Position', get(0,'Screensize'));
        sgtitle([animal ' ' region ' ' groupnames{grp} ' ITPC, black = FDR q < ' num2str(alpha)], 'Interpreter', 'none')
        for r = 1:nrow
            for sh = 1:nshank
                ch = find(chansuse == probemap(r,sh));
                if isempty(ch); continue; end
                subplot(nrow, nshank, (r-1)*nshank + sh)
                contourf(times_ms, frex, double(squeeze(itpc{grp}(ch,:,:))), 25, 'LineColor', 'none')
                if any(itpc_sig{grp}(ch,:,:), 'all') %contour of an all-false map only throws a warning
                    hold on
                    contour(times_ms, frex, double(squeeze(itpc_sig{grp}(ch,:,:))), [0.5 0.5], 'k', 'LineWidth', 1)
                    hold off
                end
                set(gca, 'YDir', 'normal', 'XLim', plot_time_ms, 'YLim', plot_freq_hz, 'CLim', clims)
                xline(0, 'r', 'LineWidth', 1)
                title(['CH' num2str(probemap(r,sh))])
            end
        end
    end
end

%% 6 CHOOSE TIME-FREQUENCY WINDOWS (after looking at sections 4/5)
%one entry per window; add more by extending the cell lists
windows = struct('name', {'win1'}, ...
                 'freq_hz', {[35 45]}, ...
                 'time_ms', {[40 75]});
nwin = numel(windows);
for w = 1:nwin
    windows(w).fidx = find(frex >= windows(w).freq_hz(1) & frex <= windows(w).freq_hz(2));
    windows(w).tidx = find(times_ms >= windows(w).time_ms(1) & times_ms <= windows(w).time_ms(2));
    if isempty(windows(w).fidx) || isempty(windows(w).tidx)
        error('Window %s falls outside frex/time_ms', windows(w).name);
    end
end

%check: windows drawn on the reference group's ITPC
plot_time_ms = [-100 300];
plot_freq_hz = [frex(1) frex(end)];
clims = [0 0.5];
probemap = ProbeInfo.ProbeMaps{prb};
[nrow, nshank] = size(probemap);
figure('Position', get(0,'Screensize'));
sgtitle([animal ' ' region ' ' groupnames{ref_group} ' ITPC with windows'], 'Interpreter', 'none')
for r = 1:nrow
    for sh = 1:nshank
        ch = find(chansuse == probemap(r,sh));
        if isempty(ch); continue; end
        subplot(nrow, nshank, (r-1)*nshank + sh)
        contourf(times_ms, frex, double(squeeze(itpc{ref_group}(ch,:,:))), 25, 'LineColor', 'none')
        set(gca, 'YDir', 'normal', 'XLim', plot_time_ms, 'YLim', plot_freq_hz, 'CLim', clims)
        for w = 1:nwin
            rectangle('Position', [windows(w).time_ms(1), windows(w).freq_hz(1), diff(windows(w).time_ms), diff(windows(w).freq_hz)], ...
                'EdgeColor', 'w', 'LineWidth', 1)
        end
        title(['CH' num2str(probemap(r,sh))])
    end
end

%% 7 PER-TRIAL PHASES INSIDE THE WINDOWS (indexed from phases, nothing recomputed)
ntot = height(trial_table);
win_u = cell(1,nwin); %complex single, nch x nfw x ntw x ntot (trial dim = trial_table rows)
ref_sum = cell(1,nwin); %sum of reference-group unit vectors, nch x nfw x ntw
for w = 1:nwin
    fi = windows(w).fidx;
    ti = windows(w).tidx;
    win_u{w} = complex(zeros(nch, numel(fi), numel(ti), ntot, 'single'));
    for grp = 1:ngrp
        rows_tt = find(trial_table.group == grp);
        win_u{w}(:,:,:,rows_tt) = exp(1i*phases{grp}(:, fi, ti, trial_table.row_in_group(rows_tt)));
    end
    ref_sum{w} = sum(win_u{w}(:,:,:,trial_table.group == ref_group), 4);
end

%% 8 WINDOW ITPC -> ONE VALUE PER CHANNEL
itpc_win = cell(1,nwin); %nch x ngrp
rayleighZ = cell(1,nwin);
for w = 1:nwin
    for grp = 1:ngrp
        itpc_win{w}(:,grp) = mean(reshape(itpc{grp}(:, windows(w).fidx, windows(w).tidx), nch, []), 2);
        rayleighZ{w}(:,grp) = ntr(grp) * itpc_win{w}(:,grp).^2;
    end
end

%optional permutation test: window ITPC difference between compare_groups, per channel
perm_p = cell(1,nwin);
perm_q = cell(1,nwin);
if do_signif == 1 && ngrp >= 2
    g1 = compare_groups(1); g2 = compare_groups(2);
    rows12 = find(ismember(trial_table.group, [g1 g2]));
    lab = trial_table.group(rows12) == g1;
    rng(1);
    for w = 1:nwin
        disp(['Permutation test: window ' windows(w).name])
        U = reshape(win_u{w}(:,:,:,rows12), nch, [], numel(rows12)); %nch x pixels x trials
        obs = mean(abs(mean(U(:,:,lab),3)),2) - mean(abs(mean(U(:,:,~lab),3)),2);
        nullD = zeros(nch, nperm);
        for k = 1:nperm
            pl = lab(randperm(numel(lab)));
            nullD(:,k) = mean(abs(mean(U(:,:,pl),3)),2) - mean(abs(mean(U(:,:,~pl),3)),2);
        end
        perm_p{w} = (sum(abs(nullD) >= abs(obs), 2) + 1) / (nperm + 1);
        [ps, si] = sort(perm_p{w});
        m = numel(ps);
        q = ps .* m ./ (1:m)';
        q = flipud(cummin(flipud(q)));
        perm_q{w} = zeros(m,1);
        perm_q{w}(si) = min(q, 1);
    end
    clear U nullD
end

%long table for R: one row per window x group x channel
wt = {};
for w = 1:nwin
    for grp = 1:ngrp
        pp = nan(nch,1); qq = nan(nch,1);
        if ~isempty(perm_p{w}); pp = perm_p{w}; qq = perm_q{w}; end
        wt{end+1} = table(repmat(string(animal),nch,1), repmat(string(region),nch,1), repmat(string(windows(w).name),nch,1), ...
            repmat(windows(w).freq_hz(1),nch,1), repmat(windows(w).freq_hz(2),nch,1), ...
            repmat(windows(w).time_ms(1),nch,1), repmat(windows(w).time_ms(2),nch,1), ...
            repmat(grp,nch,1), repmat(string(groupnames{grp}),nch,1), chansuse(:), ch_shank, ch_row, ...
            double(itpc_win{w}(:,grp)), double(rayleighZ{w}(:,grp)), repmat(ntr(grp),nch,1), pp, qq, ...
            'VariableNames', {'animal','region','window','freq_lo','freq_hi','time_lo','time_hi', ...
            'group','groupname','channel','shank','depth_row','itpc','rayleighZ','ntrials','perm_p','perm_q'}); %#ok<SAGROW>
    end
end
win_table = vertcat(wt{:});

%% 8b PLOT: window ITPC on the probe grid, per-shank depth profiles, shank means
probemap = ProbeInfo.ProbeMaps{prb};
[nrow, nshank] = size(probemap);
pr = zeros(nch,1); ps_ = zeros(nch,1);
for ch = 1:nch
    [pr(ch), ps_(ch)] = find(probemap == chansuse(ch));
end
for w = 1:nwin
    vals = itpc_win{w};
    titles = groupnames;
    if size(vals,2) >= 2
        vals(:,end+1) = vals(:,compare_groups(1)) - vals(:,compare_groups(2));
        titles{end+1} = [groupnames{compare_groups(1)} ' - ' groupnames{compare_groups(2)}];
    end
    figure('Position', [100 100 300*size(vals,2) 500]);
    sgtitle([animal ' ' region ' window ' windows(w).name ' ITPC'], 'Interpreter', 'none')
    for k = 1:size(vals,2)
        mat = nan(nrow, nshank);
        mat(sub2ind([nrow nshank], pr, ps_)) = vals(:,k);
        subplot(1, size(vals,2), k)
        imagesc(mat, 'AlphaData', ~isnan(mat))
        colorbar
        title(titles{k}, 'Interpreter', 'none')
        xlabel('Shank'); ylabel('Depth row')
        for ch = 1:nch
            lbl = num2str(chansuse(ch));
            if k == size(vals,2) && ~isempty(perm_q{w}) && perm_q{w}(ch) < alpha
                lbl = [lbl '*']; %significant group difference
            end
            text(ps_(ch), pr(ch), lbl, 'HorizontalAlignment', 'center', 'FontSize', 7)
        end
    end

    %depth profile per shank, one subplot per group
    figure('Position', [100 100 350*ngrp 400]);
    sgtitle([animal ' ' region ' window ' windows(w).name ' ITPC by depth'], 'Interpreter', 'none')
    for grp = 1:ngrp
        subplot(1, ngrp, grp)
        hold on
        for sh = 1:nshank
            sel = find(ps_ == sh);
            [rr, o] = sort(pr(sel));
            plot(rr, itpc_win{w}(sel(o),grp), '-o', 'DisplayName', ['Shank ' num2str(sh)])
        end
        hold off
        legend('Location', 'best')
        title(groupnames{grp}, 'Interpreter', 'none')
        xlabel('Depth row'); ylabel('Window ITPC')
    end

    %shank means
    shank_mean = nan(nshank, ngrp);
    for sh = 1:nshank
        shank_mean(sh,:) = mean(itpc_win{w}(ps_ == sh,:), 1);
    end
    figure();
    bar(shank_mean)
    legend(groupnames, 'Interpreter', 'none')
    xlabel('Shank'); ylabel('Mean window ITPC')
    title([animal ' ' region ' window ' windows(w).name ' shank means'], 'Interpreter', 'none')
end

%% 9 PER-TRIAL PHASE-RESET SCORE (metric PROVISIONAL - to be discussed)
reset_method = 'alignment'; %'alignment' | 'consistency' | 'circdist'
isref = trial_table.group == ref_group;
reset_d = cell(1,nwin); %trial phase - reference angle, radians (-pi, pi]
reset_score = cell(1,nwin); %nch x ntot
reset_offset = cell(1,nwin); %nch x ntot, circular mean of reset_d over the window
for w = 1:nwin
    refv = repmat(ref_sum{w}, 1, 1, 1, ntot);
    refv(:,:,:,isref) = refv(:,:,:,isref) - win_u{w}(:,:,:,isref); %leave-one-out for reference trials
    reset_d{w} = angle(win_u{w} .* conj(refv));
    D = reshape(reset_d{w}, nch, [], ntot);
    z = mean(exp(1i*D), 2);
    reset_offset{w} = reshape(angle(z), nch, ntot);
    switch reset_method
        case 'alignment' %mean cos: -1..1, 1 = locked to the reference phase
            reset_score{w} = reshape(mean(cos(D), 2), nch, ntot);
        case 'consistency' %resultant length: 0..1, catches a reset at a constant phase offset
            reset_score{w} = reshape(abs(z), nch, ntot);
        case 'circdist' %mean absolute phase distance: 0..pi, 0 = locked to the reference phase
            reset_score{w} = reshape(mean(abs(D), 2), nch, ntot);
    end
end
clear refv D z

%long table for R: one row per window x trial x channel
rt = {};
for w = 1:nwin
    base = trial_table(repelem((1:ntot)', nch), :);
    base.window = repmat(string(windows(w).name), height(base), 1);
    base.channel = repmat(chansuse(:), ntot, 1);
    base.shank = repmat(ch_shank, ntot, 1);
    base.depth_row = repmat(ch_row, ntot, 1);
    base.reset_score = double(reshape(reset_score{w}, [], 1));
    base.reset_offset = double(reshape(reset_offset{w}, [], 1));
    base.reset_method = repmat(string(reset_method), height(base), 1);
    rt{end+1} = base; %#ok<SAGROW>
end
trial_score_table = vertcat(rt{:});
trial_score_table = [table(repmat(string(animal),height(trial_score_table),1), repmat(string(region),height(trial_score_table),1), ...
    'VariableNames', {'animal','region'}), trial_score_table];

%% 9b PLOT: reset scores
chcheck = chansuse(1); %raw channel ID for the distribution plot
probemap = ProbeInfo.ProbeMaps{prb};
pr = zeros(nch,1); ps_ = zeros(nch,1);
for ch = 1:nch
    [pr(ch), ps_(ch)] = find(probemap == chansuse(ch));
end
[~, ord] = sortrows([ps_ pr]); %channels ordered by shank, then depth
grp_of_trial = trial_table.group;
for w = 1:nwin
    figure('Position', [100 100 1200 500]);
    sgtitle([animal ' ' region ' window ' windows(w).name ' reset score (' reset_method ')'], 'Interpreter', 'none')

    %trial x channel matrix: does a reset span the probe?
    subplot(1,2,1)
    imagesc(reset_score{w}(ord,:))
    colorbar
    hold on
    bounds = find(diff(grp_of_trial)) + 0.5;
    for b = bounds'
        xline(b, 'w', 'LineWidth', 1.5);
    end
    shb = find(diff(ps_(ord))) + 0.5;
    for b = shb'
        yline(b, 'w', 'LineWidth', 1);
    end
    hold off
    set(gca, 'YTick', 1:nch, 'YTickLabel', chansuse(ord), 'FontSize', 6)
    xlabel('Trial (grouped by group)'); ylabel('Channel (shank, depth)')

    %score distribution per group for one channel
    subplot(1,2,2)
    ch = find(chansuse == chcheck);
    hold on
    for grp = 1:ngrp
        histogram(reset_score{w}(ch, grp_of_trial == grp), 20, 'DisplayName', groupnames{grp})
    end
    hold off
    legend('Interpreter', 'none')
    xlabel('Reset score'); ylabel('Trials')
    title(['CH' num2str(chcheck)])
end

%% 10 ISPC: channel x channel phase clustering inside each window, per group
ispc = cell(1,nwin);
for w = 1:nwin
    for grp = 1:ngrp
        rows_g = trial_table.group == grp;
        U = reshape(win_u{w}(:,:,:,rows_g), nch, [], sum(rows_g)); %nch x pixels x trials
        acc = zeros(nch);
        for p = 1:size(U,2)
            X = double(reshape(U(:,p,:), nch, [])); %nch x trials
            acc = acc + abs(X*X') / size(X,2); %|mean over trials of exp(1i*(phiA - phiB))|
        end
        ispc{w}{grp} = acc / size(U,2);
    end
end
clear U X

%% 10b PLOT: ISPC heatmaps (ordered by shank, then depth)
probemap = ProbeInfo.ProbeMaps{prb};
pr = zeros(nch,1); ps_ = zeros(nch,1);
for ch = 1:nch
    [pr(ch), ps_(ch)] = find(probemap == chansuse(ch));
end
[~, ord] = sortrows([ps_ pr]);
shb = find(diff(ps_(ord))) + 0.5;
for w = 1:nwin
    mats = ispc{w};
    titles = groupnames;
    clim_list = repmat({[0 1]}, 1, numel(mats));
    if numel(mats) >= 2
        mats{end+1} = ispc{w}{compare_groups(1)} - ispc{w}{compare_groups(2)};
        titles{end+1} = [groupnames{compare_groups(1)} ' - ' groupnames{compare_groups(2)}];
        clim_list{end+1} = [-1 1] * max([abs(mats{end}(:)); eps]);
    end
    figure('Position', [100 100 450*numel(mats) 420]);
    sgtitle([animal ' ' region ' window ' windows(w).name ' ISPC'], 'Interpreter', 'none')
    for k = 1:numel(mats)
        subplot(1, numel(mats), k)
        imagesc(mats{k}(ord,ord), clim_list{k})
        axis square
        colorbar
        hold on
        for b = shb'
            xline(b, 'w', 'LineWidth', 1);
            yline(b, 'w', 'LineWidth', 1);
        end
        hold off
        set(gca, 'XTick', 1:nch, 'XTickLabel', chansuse(ord), 'YTick', 1:nch, 'YTickLabel', chansuse(ord), 'FontSize', 5)
        title(titles{k}, 'Interpreter', 'none')
    end
end

%% 11 OUTPUT STRUCT + RUN REPORT
itpc_info = struct();
itpc_info.ITPC = itpc;
itpc_info.prefAngles = prefAngle;
itpc_info.dimfeatures = {'chans','freqs','time'};
itpc_info.stimtime = find(times_ms == 0);
itpc_info.frex = frex;
itpc_info.times_ms = times_ms;
itpc_info.range_cycles = range_cycles;
itpc_info.windows = windows;
itpc_info.ch_index = chansuse;
itpc_info.ch_shank = ch_shank;
itpc_info.ch_row = ch_row;
itpc_info.itpc_win = itpc_win;
itpc_info.rayleighZ = rayleighZ;
itpc_info.itpc_p = itpc_p;
itpc_info.itpc_sig = itpc_sig;
itpc_info.perm_p = perm_p;
itpc_info.perm_q = perm_q;
itpc_info.compare_groups = compare_groups;
itpc_info.ref_group = ref_group;
itpc_info.trial_table = trial_table;
itpc_info.win_u = win_u;
itpc_info.reset_d = reset_d;
itpc_info.reset_score = reset_score;
itpc_info.reset_offset = reset_offset;
itpc_info.reset_method = reset_method;
itpc_info.ispc = ispc;
itpc_info.condition_directory = condi_dir;
itpc_info.all_conditions = condis;
itpc_info.condition_grouping = cellfun(@(g) condis(g), groupassign, 'UniformOutput', false);
itpc_info.condition_grouping_index = groupassign;
itpc_info.groupnames = groupnames;
itpc_info.region = region;
itpc_info.null_check = null_check;
itpc_info.report = report;

fprintf('\n===== ITPC_PhaseReset_singleanimal run report: %s %s =====\n', animal, region);
fprintf('Channels: %d used of %d on probe (Ch_Remove: %s)\n', report.nchans_used, report.nchans_probe, mat2str(report.ch_removed));
fprintf('Trials per group: %s\n', strjoin(compose('%s = %d', string(groupnames(:)), ntr(:)), ', '));
fprintf('Stored phases: %.1f GB\n', report.phase_gb);
disp(report.files)
if isfield(report, 'sig_channels')
    fprintf('Channels with any FDR-significant post-stim ITPC pixel: %s\n', strjoin(compose('%s = %d', string(groupnames(:)), report.sig_channels(:)), ', '));
end
for w = 1:nwin
    if ~isempty(perm_q{w})
        fprintf('Window %s: %d channels with q < %g (%s vs %s)\n', windows(w).name, sum(perm_q{w} < alpha), alpha, ...
            groupnames{compare_groups(1)}, groupnames{compare_groups(2)});
    end
end
fprintf('Reset score method: %s (PROVISIONAL)\n', reset_method);
for k = 1:numel(report.notes)
    fprintf('NOTE: %s\n', report.notes{k});
end

%% 12 SAVE (itpc_info .mat + CSVs for R; phases are not saved)
folder = fullfile(save_root, animal);
if ~isfolder(folder)
    mkdir(folder)
end
base = [animal '_' region '_' strjoin(groupnames, '_vs_')];
if null_check == 1
    base = [base '_NULLCHECK'];
end
savefile = fullfile(folder, [base '_itpc.mat']);
dosave = true;
if isfile(savefile)
    answer = questdlg(['Overwrite ' base '_itpc.mat and its CSVs?'], 'File exists', 'Overwrite', 'Cancel', 'Cancel');
    dosave = strcmp(answer, 'Overwrite');
end
if dosave
    save(savefile, 'itpc_info', '-v7.3')
    writetable(win_table, fullfile(folder, [base '_windowITPC.csv']))
    writetable(trial_score_table, fullfile(folder, [base '_trialReset.csv']))
    disp(['Saved to ' folder])
else
    disp('Not saved')
end
