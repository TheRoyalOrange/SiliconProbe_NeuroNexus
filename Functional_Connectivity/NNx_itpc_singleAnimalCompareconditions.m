% itpc_singleAnimalCompareconditions.m
%
% Description: For one animal, loads the saved ITPC results of two
%   conditions (one file per condition, from the save line at the bottom of
%   ITPC_PhaseReset_singleanimal.m), takes the difference condition 1 minus
%   condition 2 of the measure you pick (ITPC, Rayleigh Z or 2nd-harmonic
%   ITPC), and plots that difference as a time-frequency map per channel in
%   the probe layout, one figure per probe. With sig_only = 1 the
%   difference is only shown where the reference condition (ref_cond) has
%   significant ITPC (Rayleigh p < sig_p per point, uncorrected, threshold
%   from its trial count), the rest blanked. Optionally draws the ERPs of
%   both conditions on top (right-hand y axis); then it asks for an ERP
%   band (inputdlg, empty = raw ERP), bandpasses both ERPs, redraws, and
%   repeats until you're happy (as in section 11 of the ITPC script).
%
% Inputs:
%   ProbeInfo (struct) - loaded from <probeinfo_dir>\<animal>-ProbeInfo.mat
%     .Areas (cell of char, 1 x probenum) - region per probe; regions picks
%       from these ({} = all)
%     .ProbeIds (cell, 1 x probenum) - raw channel IDs per probe
%     .Ch_Remove (cell, 1 x probenum) - raw channel IDs to exclude
%     .ProbeMaps (cell, 1 x probenum) - 2D [depth rows x shanks] matrix of
%       raw channel IDs; columns = shanks
%   one ITPC results .mat per condition (picked with uipickfiles), saved by
%     ITPC_PhaseReset_singleanimal.m, containing:
%     itpc_lfp (single, nch x nfrex x ntime) - ITPC, 0-1
%     itpc2_lfp (single, nch x nfrex x ntime) - 2nd-harmonic ITPC, 0-1
%     rayleighZ_lfp (single, nch x nfrex x ntime) - ntrials * ITPC^2
%     erp_raw (double, nseg x nch) - trial-average LFP, unfiltered, uV
%   INFERRED (the files don't store these, so they're rebuilt here and must
%   match the ITPC run that made the files):
%     - channel dim k = chansuse(k), built from ProbeInfo with the same
%       regions and the same Ch_Remove as that run (probes stacked in order)
%     - frequency dim = frex, time dim = time_ms(1):time_ms(2) (1 ms steps)
%     - erp_raw time = (time_ms(1)-pad_ms):(time_ms(2)+pad_ms), ms re: stim
%     Sizes are checked against these and the script stops if they don't fit.
%   User parameters - see the Setup section.
%
% Outputs (workspace):
%   res (cell, 1 x 2) - the loaded files, one struct per condition with the
%     four variables above
%   itpc_diff (single, nch x nfrex x ntime) - condition 1 minus condition 2
%     of the measure picked by diff_measure (diff_name says which), unmasked
%   sigmask (logical, nch x nfrex x ntime) - true where ref_cond's Rayleigh
%     Z >= rayZ_sig (p < sig_p)
%   itpc_diff_show (single, same size) - itpc_diff with NaN outside sigmask
%     when sig_only = 1 (what gets plotted)
%   erp_plot (cell, 1 x 2) - (double, nch x nseg) ERP per condition, uV,
%     bandpassed to erp_band if one is set
%   erp_band (double, 1 x 2 or empty) - last band picked, Hz ([] = raw)
%   chansuse, ch_probe, ch_row, ch_shank (double, nch x 1) - raw channel ID,
%     probe index, ProbeMaps row and shank of each channel
%   frex, times_ms, seg_ms - axes of the frequency, TF time and ERP time dims
%   report (struct) - files used, estimated trial count per condition
%     (rayleighZ / itpc^2), notes; shown in the last section
%
% Dependencies: ITPC_PhaseReset_singleanimal.m (results files, with its
%   save line), ProbeInfo from NNx_BaseAnalysis.m,
%   uipickfiles.m, Signal Processing Toolbox (bandpass, only if erp_band
%   is set).

%% Setup: Specify Animal Info, Conditions, and Plot Settings

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Who is the data?
animal = '20261011-p7d';

%which probes were in the ITPC run? MUST be the same regions (and the same ProbeInfo Ch_Remove) as when the ITPC files were saved
regions = {}; %{} = every probe in the recording

%which two conditions are we comparing? (you get one file-picker prompt each, in this order)
condnames = {'WhiskerActive', 'WhiskerActive_TTX'}; %the difference is condition 1 minus condition 2

%where are the ITPC results files?
load_directory = pwd; %the ITPC script saves into whatever folder Matlab was in

%what are we taking the difference of?
diff_measure = 1;
    % 1 = ITPC (biased upward for small trial counts, so careful if the conditions have very different n)
    % 2 = Rayleigh Z (n * ITPC^2, ~1 with no phase locking whatever the trial count)
    % 3 = 2nd-harmonic ITPC (trials split into two groups pi apart)

%which condition is the reference? (index into condnames) only points where it has significant ITPC get shown
ref_cond = 1;
sig_only = 0;
    % 0 = show the whole difference
    % 1 = blank (NaN) everything outside the reference condition's Rayleigh p < sig_p (outlined in grey)
sig_p = 0.05; %per point, uncorrected (same as the section 10b plot in the ITPC script)

%%%ITPC run settings (copy these from section 1 of ITPC_PhaseReset_singleanimal.m, the files don't store them)%%%
frex = 1:1:80; %Hz
time_ms = [-1000 1000]; %ms re: stim
pad_ms = 1500; %ms

%%%Plotting details%%%

%what time and frequency range to show?
itpcplot_ms = [-200 500]; %ms re: stim
itpcplot_hz = [1 80]; %Hz
diff_clim = []; %colour limits for the difference, e.g. [-0.2 0.2]; [] = symmetric around 0, from each probe's biggest difference

%want the ERPs of both conditions drawn on top?
plot_erp = 0;
    % 0 = difference only
    % 1 = ERPs on a right-hand y axis (condition 1 solid black, condition 2 dashed black)
erp_band = []; %Hz, the band offered first in the band picker (you get asked anyway), e.g. [4 8]; [] = raw ERP

%save the figures?
plotsave = 0;
    % 0 = no
    % 1 = yes, into save_directory\Figures_ITPC\animal
save_directory = 'E:\Roy\Silcon Probe Data Followup Analyses\ITPC';

probeinfo_dir = 'E:\Roy\Processed Silicon Probe Data\ProbeInfo';

%%%nothing else to fill in%%%

%% Data Extraction: ProbeInfo and the channels in the ITPC files
report = struct();
report.notes = {};

load(fullfile(probeinfo_dir, [animal '-ProbeInfo.mat']), 'ProbeInfo');

%which probes? (same picking as the ITPC script)
if isempty(regions) == true
    prbs = 1:numel(ProbeInfo.Areas); %nothing picked = every probe
else
    prbs = find(ismember(ProbeInfo.Areas, regions)); %just the ones you named
end

%channels: raw IDs on each probe minus its Ch_Remove, stacked probe after probe (same as the ITPC script).
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
end
nch = numel(chansuse); %all probes together

%time axes of the saved arrays
times_ms = time_ms(1):time_ms(2); %ms re: stim, time dim of the ITPC arrays
seg_ms = (time_ms(1)-pad_ms):(time_ms(2)+pad_ms); %ms re: stim, time dim of erp_raw

%% Data Extraction: load the ITPC results of each condition
res = cell(1,2); %one struct per condition
for c = 1:2
    pickedfile = uipickfiles('FilterSpec', load_directory, 'Prompt', ['Choose the ' condnames{c} ' ITPC results file (just one)'], 'NumFiles', 1);
    if ~iscell(pickedfile) || isempty(pickedfile)
        error('No file chosen for condition %s', condnames{c});
    end
    disp(['Loading: ' condnames{c}])
    res{c} = load(pickedfile{1}, 'itpc_lfp', 'itpc2_lfp', 'rayleighZ_lfp', 'erp_raw');
    report.files{c} = pickedfile{1};

    %do the sizes fit the channels and axes rebuilt above? if not, the setup doesn't match the ITPC run
    if isequal(size(res{c}.itpc_lfp), [nch numel(frex) numel(times_ms)]) == false
        error('%s: itpc_lfp is %s but the setup gives [%d %d %d] (channels x freqs x time). Check regions, Ch_Remove, frex and time_ms', ...
            condnames{c}, mat2str(size(res{c}.itpc_lfp)), nch, numel(frex), numel(times_ms));
    end
    if isequal(size(res{c}.erp_raw), [numel(seg_ms) nch]) == false
        error('%s: erp_raw is %s but the setup gives [%d %d] (time x channels). Check time_ms and pad_ms', ...
            condnames{c}, mat2str(size(res{c}.erp_raw)), numel(seg_ms), nch);
    end

    %trial count isn't saved, but rayleighZ = n * ITPC^2, so we can get it back
    ok = res{c}.itpc_lfp > 0.01; %skip ~0 ITPC points (dividing by ~0)
    report.ntrials(c) = round(median(res{c}.rayleighZ_lfp(ok) ./ res{c}.itpc_lfp(ok).^2));
end
report.condnames = condnames;

%% Data PreProcessing: condition 1 minus condition 2
diff_fields = {'itpc_lfp', 'rayleighZ_lfp', 'itpc2_lfp'}; %what diff_measure picks
diff_labels = {'ITPC', 'Rayleigh Z', '2nd-harmonic ITPC'};
diff_name = diff_labels{diff_measure};
itpc_diff = res{1}.(diff_fields{diff_measure}) - res{2}.(diff_fields{diff_measure}); % [channels x freqs x time]

if diff_measure == 1 && report.ntrials(1) ~= report.ntrials(2)
    %chance ITPC is ~sqrt(pi/(4n)), so different trial counts shift the ITPC difference on their own
    report.notes{end+1} = sprintf('Trial counts differ (%d vs %d): chance ITPC is %.3f vs %.3f, so part of the ITPC difference is just n. Rayleigh Z (diff_measure = 2) is n-corrected.', ...
        report.ntrials(1), report.ntrials(2), sqrt(pi/(4*report.ntrials(1))), sqrt(pi/(4*report.ntrials(2))));
end

%where is the reference condition significant? Rayleigh Z at p = sig_p for its trial count:
%the Zar formula from the ITPC script, p = exp(sqrt(1 + 4n + 4(n^2 - Rn^2)) - (1 + 2n)), solved for Z = Rn^2/n
ntr_ref = report.ntrials(ref_cond);
L = log(sig_p) + 1 + 2*ntr_ref;
rayZ_sig = (1 + 4*ntr_ref + 4*ntr_ref^2 - L^2) / (4*ntr_ref); %Rayleigh Z threshold
sigmask = res{ref_cond}.rayleighZ_lfp >= rayZ_sig; % [channels x freqs x time], true = reference significant
report.sig_frac = mean(sigmask(:)); %fraction of points that survive (all channels, freqs, times)

itpc_diff_show = itpc_diff; %what gets plotted; itpc_diff itself stays whole
if sig_only == 1
    itpc_diff_show(sigmask == 0) = NaN; %blank everything the reference isn't significant at
end

%% Data Visualization: difference in the probe layout (+ ERPs if you want)
%with plot_erp = 1 you get asked for an ERP band (empty = raw ERP), and it repeats until you're happy, like section 11 of the ITPC script
disp('Plotting: ITPC difference')
tpi = find(times_ms >= itpcplot_ms(1) & times_ms <= itpcplot_ms(2)); %time points to show
fpi = find(frex >= itpcplot_hz(1) & frex <= itpcplot_hz(2)); %frequencies to show
eidx = find(seg_ms >= itpcplot_ms(1) & seg_ms <= itpcplot_ms(2)); %ERP points that line up with the TF time axis

%blue-white-red colormap: blue = condition 2 higher, red = condition 1 higher
ncol = 128;
bwr = [linspace(0,1,ncol)' linspace(0,1,ncol)' ones(ncol,1); ...
       ones(ncol,1) linspace(1,0,ncol)' linspace(1,0,ncol)'];

%create figure folder
if plotsave == 1 && isfolder([save_directory '\Figures_ITPC\' animal]) == false
    mkdir([save_directory '\Figures_ITPC\' animal])
end

ishsokay = 0;
while ishsokay == 0
    %which band for the ERPs? (only asked when they're being drawn)
    if plot_erp == 1
        dlgtitle = 'Band ERP';
        promt = {'Which frequency band for the ERPs? ([low high] in Hz, leave empty for the raw ERP)'};
        fieldsize = [1 50];
        definput = {num2str(erp_band)};
        opts.Resize = 'on';
        opts.WindowStyle = 'normal';
        bandpick = inputdlg(promt,dlgtitle,fieldsize,definput,opts);
        erp_band = str2num(bandpick{1}); %Hz, [] = raw
    end

    %ERPs to draw on top, [channels x segment]
    %the ERP is already a trial average and the bandpass is linear and zero-phase, so filtering the average = averaging filtered trials
    erp_plot = cell(1,2);
    for c = 1:2
        if isempty(erp_band) == true
            erp_plot{c} = res{c}.erp_raw'; %raw trial-average, uV
        else
            erp_plot{c} = bandpass(res{c}.erp_raw, erp_band, 1000, 'ImpulseResponse', 'iir')'; %zero-phase bandpass, uV (bandpass works down columns = time)
        end
    end
    if isempty(erp_band) == true
        erp_label = 'raw ERP';
    else
        erp_label = sprintf('%g-%g Hz ERP', erp_band);
    end

    for prb = prbs
        probemap = ProbeInfo.ProbeMaps{prb}; % [depth rows x shanks]
        [nrow, nshank] = size(probemap);
        onprobe = ch_probe == prb;

        %same colour scale for the whole probe, centred on 0
        if isempty(diff_clim) == true
            dlim = max(abs(itpc_diff_show(onprobe, fpi, tpi)), [], 'all'); %max skips the blanked NaNs
            if isnan(dlim) == true || dlim == 0
                dlim = 1; %nothing survived on this probe, any limits will do
            end
            clim_p = [-dlim dlim];
        else
            clim_p = diff_clim;
        end
        %same ERP scale for the whole probe, both conditions
        elim = max([max(abs(erp_plot{1}(onprobe, eidx)), [], 'all') max(abs(erp_plot{2}(onprobe, eidx)), [], 'all')]);

        figure();
        set(gcf, 'Position', get(0, 'Screensize'));
        ttl = sprintf('%s %s-Probe %s: %s minus %s', animal, ProbeInfo.Areas{prb}, diff_name, condnames{1}, condnames{2});
        if sig_only == 1
            ttl = [ttl sprintf(', only where %s ITPC is significant (Rayleigh p < %g, grey outline)', condnames{ref_cond}, sig_p)];
        end
        if plot_erp == 1
            ttl = [ttl sprintf('. ERPs (right axis, uV, %s): %s solid, %s dashed', erp_label, condnames{1}, condnames{2})];
        end
        sgtitle(ttl, 'Interpreter', 'none')

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
                tf = reshape(itpc_diff_show(k, fpi, tpi), numel(fpi), numel(tpi)); % [frequencies x time]
                sm = reshape(sigmask(k, fpi, tpi), numel(fpi), numel(tpi)); %reference significant, same shape

                if plot_erp == 1
                    yyaxis left %difference on the left axis
                end
                if any(isnan(tf(:)) == false) %only fill if something survived (contour complains about all-NaN panels)
                    contour(times_ms(tpi), frex(fpi), tf, 50, 'Fill', 'on')
                end
                if sig_only == 1 && any(sm(:)) && all(sm(:)) == false
                    hold on
                    contour(times_ms(tpi), frex(fpi), double(sm), [0.5 0.5], 'LineColor', [0.4 0.4 0.4], 'LineWidth', 1); %edge of the significant region
                    hold off
                end
                set(ax, 'YDir', 'normal', 'CLim', clim_p, 'Box', 'on')
                colormap(ax, bwr)
                ylim([frex(fpi(1)) frex(fpi(end))])

                if plot_erp == 1
                    yyaxis right %ERPs on the right axis
                    hold on
                    plot(seg_ms(eidx), erp_plot{1}(k, eidx), 'k', 'LineWidth', 1)
                    plot(seg_ms(eidx), erp_plot{2}(k, eidx), 'k', 'LineWidth', 1, 'LineStyle','--')
                    hold off
                    if elim > 0
                        ylim([-elim elim])
                    end
                    ax.YAxis(1).Color = 'k'; %both axes black instead of yyaxis' default blue/orange
                    ax.YAxis(2).Color = 'k';
                end
                xline(0, 'LineWidth', 1, 'Color', 'g'); %stim onset (green, red would get lost in the colormap)
                xlim(itpcplot_ms)

                title(['CH' num2str(chid)])
                if row < nrow
                    set(ax, 'XTickLabel', [])
                end
                if row == nrow && shank == 1
                    xlabel('Time (ms)')
                    if plot_erp == 1
                        yyaxis left
                        ylabel('Frequency (Hz)')
                        yyaxis right
                        ylabel('ERP (uV)')
                    else
                        ylabel('Frequency (Hz)')
                    end
                end
                lastax = ax;
            end
        end
        if isempty(lastax) == false
            axpos = lastax.Position;
            cb = colorbar(lastax, 'eastoutside');
            cb.Label.String = [diff_name ' difference'];
            lastax.Position = axpos; %keep the panel the same size as the others
        end

        if plotsave == 1
            fname = sprintf([animal '-' ProbeInfo.Areas{prb} '_Probe' '-' strrep(diff_name, ' ', '') 'diff_' condnames{1} '_minus_' condnames{2}]);
            if sig_only == 1
                fname = [fname '_sig' condnames{ref_cond}]; %masked version gets its own file
            end
            if plot_erp == 1
                fname = [fname '_' strrep(erp_label, ' ', '')]; %band in the name, so each band you try gets its own file
            end
            savefig(gcf, [save_directory '\Figures_ITPC\' animal '\' fname], 'compact')
        end
    end

    %are you happy with that? if yes, continue. If no, pick another band
    if plot_erp == 1
        quest = 'Happy with that ERP band?';
        btn1 = 'Yes, happy as a clam';
        btn2 = 'No, lemme pick again';
        answer = questdlg(quest, 'Band ERP', btn1, btn2, btn1);
        switch answer
            case btn1
                ishsokay = 1;
            case btn2
                ishsokay = 0;
        end
    else
        ishsokay = 1; %no ERPs, nothing to pick, so one pass is it
    end
end

%% run report
disp('===== itpc_singleAnimalCompareconditions run report =====')
for c = 1:2
    disp([condnames{c} ': ~' num2str(report.ntrials(c)) ' trials, file ' report.files{c}])
end
disp(sprintf('Reference %s significant (Rayleigh Z >= %.2f, p < %g, n ~ %d) at %.1f%% of channel x freq x time points', ...
    condnames{ref_cond}, rayZ_sig, sig_p, ntr_ref, 100*report.sig_frac)) %#ok<DSPSP>
for k = 1:numel(report.notes)
    disp(['NOTE: ' report.notes{k}])
end
