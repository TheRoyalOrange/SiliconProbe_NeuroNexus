%% OscillationEvents_StimPrePost.m
%
% Description: Runs oshkosh_detect_funct on every recording you pick for one
%   animal, sorts the probe events (prbEvent) by whether a stimulus was
%   present, and measures the participation ratio (PR) of the filtered LFP
%   across each probe's channels, for each event and in windows just before and
%   just after each stimulus. PR = (sum of eigenvalues)^2 / (sum of squared
%   eigenvalues) of the channel covariance matrix. It's about how many
%   independent dimensions the activity uses: 1 = every channel doing the same
%   thing, number of channels = every channel doing its own thing. Also opens
%   the "all events" figures (made inside oshkosh_detect_funct, one per
%   recording) and the scrollable event viewer (oscillation_plotEvents).
%   Nothing is saved to disk, everything stays in the workspace.
%
% Inputs:
%   animal (string) - animal ID, e.g. '20260821-p10'. Recording folder names
%     must contain <animal>-<condition> (oshkosh_detect_funct errors otherwise)
%   save_directory (string) - folder holding <animal>-ProbeInfo.mat
%   fq_range (double, 1x2, Hz) - bandpass used for oscillation detection, and
%     so also the band PR is measured in
%   probes (int vector) - indices into ProbeInfo.ProbeIds/ProbeMaps to analyze
%   reorder (0/1) - pre-06/2024 channel remap flag for oshkosh_detect_funct
%   nostim_conds (cell of strings) - recordings whose condition name contains
%     one of these have their stim_times and during_stim flags wiped
%   stim_reach (double, ms) - an event starting this long after a stim still
%     counts as "with stim". 3000 matches the viewer's grey shading
%   train_gap (double, ms) - TTL-high samples closer than this are one stim
%   pre_win, post_win (double, 1x2, ms) - windows relative to stim onset
%   pr_matrix (0/1) - build PR from channel covariance or correlation
%   view_which (0/1/2) - which events the scrollable viewer shows
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
%   eventTable (cell, 1 x length(probes)) - each a table, one row per prbEvent:
%     file, condition, start, stop (samples), dur_ms, peak_rms (uV),
%     during_stim (TTL overlap only), with_stim (TTL overlap or within
%     stim_reach), pr (participation ratio over the event, 1 to Nch)
%   stimSortSummary (cell, 1 x length(probes)) - table, one row per stim group
%     (no stim / with stim): n_events and medians of the eventTable measures
%   pr_pre, pr_post (cell, 1 x length(probes)) - each Nstim x 1: participation
%     ratio in that window for each stim, NaN if no event snippet covers it
%   prepostSummary (cell, 1 x length(probes)) - table, one row per stim
%     condition: n_stim, n_used (stims with BOTH windows covered), and mean
%     pr_pre, pr_post over those n_used stims
%   Figures - one "all events" figure per recording, one summary figure per
%     probe, one scrollable viewer
%
% Dependencies: NNx_BaseAnalysis.m (builds ProbeInfo),
%   NNx_editProbeInfo_ChRemove.m (Ch_Remove), oshkosh_detect_funct.m,
%   oscillation_plotEvents.m, uipickfiles, Signal Processing Toolbox.
%
% NOTE (inferred data contracts, not verified on real recordings):
%   - oscillations.stim_times holds EVERY TTL-high sample, not onsets, so
%     onsets are rebuilt here from gaps in that list.
%   - prbsnippets.event_data{p}{k} rows are in oscillations.chanals order, and
%     its columns line up with oscillations.prbEvent_ext_timestamps{p}(k,:).
%   - oshkosh_detect_funct only keeps signal around events (event +/- 1s), so a
%     stim window only gets a PR if an event's snippet covers all of it. Stims
%     with no event nearby are NaN, so pr_pre and pr_post describe stims near
%     an event, not all stims. n_used tells you how many that was.
%   - PR depends on window length, since short windows have few independent
%     samples. pr_pre vs pr_post are fair to compare when the windows match in
%     length. eventTable.pr uses each event's own length, so check dur_ms too
%     before comparing event groups.
%   - oshkosh_detect_funct normalizes the TTL to its own max, so a recording
%     with no stim can read noise as stim (seen on synthetic data). 1ms pulses
%     can also be missed, since it keeps only 1 of every 30 samples.


%% NEEDS INPUT FIRST: SET UP ANALYSIS INFO

%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY
%THERE ARE VARIABLES IN THIS SECTION THAT YOU MUST MODIFY

%Who is the data?
animal = '20260921-p8';

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
pr_matrix = 0;
    % 0 = covariance (big-amplitude channels count for more)
    % 1 = correlation (every channel counts the same)

%which events should the scrollable viewer show?
view_which = 0;
    % 0 = all events
    % 1 = only with stim
    % 2 = only without stim

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
nstim = numel(stim_onset);
fileReport = table((1:nfiles)', file_condition, n_ttl_samples, n_onsets, ...
    'VariableNames', {'file','condition','n_ttl_samples','n_onsets'}); %just take a looksee after the run

%how to get PR from a [time x channels] chunk of data
if pr_matrix == 0
    chan_matrix = @cov;
else
    chan_matrix = @corrcoef;
end
%for a symmetric matrix, sum of eigenvalues = trace and sum of squared eigenvalues = sum of all squared entries, so no eig needed
pr_calc = @(C) trace(C)^2 / sum(C(:).^2);


%% Data PreProcessing: Sort probe events by stim present or not, and get each event's PR
eventTable = cell(1,length(probes));
stimSortSummary = cell(1,length(probes));
with_stim = cell(1,length(probes));
for prb = 1:length(probes)
    ev_times = oscillations.prbEvent_timestamps{prb}; %[events x 2], start and stop samples
    ev_ext = oscillations.prbEvent_ext_timestamps{prb}; %[events x 2], same but +/- 1s, what the snippets cover
    ev_file = oscillations.prbEvent_file_index{prb};
    during_stim = oscillations.prbEvent_during_stim{prb};
    prb_chans = find(oscillations.chan_probe == prb); %snippets hold every analyzed channel, so this picks out this probe's rows
    nev = size(ev_times,1);

    with_stim{prb} = false(nev,1);
    pr = zeros(nev,1);
    for ev = 1:nev
        %does a stim overlap this event, or come just before it?
        ttl_high = oscillations.stim_times{ev_file(ev)};
        with_stim{prb}(ev) = during_stim(ev) == 1 | any(ttl_high <= ev_times(ev,1) & ttl_high >= ev_times(ev,1)-stim_reach);

        %PR over the event itself, without the extra 1s on each side
        ev_idx = (ev_times(ev,1):ev_times(ev,2)) - ev_ext(ev,1) + 1; %recording samples to snippet samples
        ev_dat = prbsnippets.event_data{prb}{ev}(prb_chans,ev_idx); %[channels x time]
        pr(ev) = pr_calc(chan_matrix(ev_dat'));
    end

    eventTable{prb} = table(ev_file, oscillations.prbEvent_condition{prb}, ev_times(:,1), ev_times(:,2), ...
        oscillations.prbEvent_durations{prb}, oscillations.prbEvent_peakNormedPower{prb}, ...
        logical(during_stim), with_stim{prb}, pr, ...
        'VariableNames', {'file','condition','start','stop','dur_ms','peak_rms','during_stim','with_stim','pr'});

    %medians for each group
    group_names = {'no stim'; 'with stim'};
    n_events = zeros(2,1);
    med_dur_ms = zeros(2,1);
    med_peak_rms = zeros(2,1);
    med_pr = zeros(2,1);
    for g = 1:2
        in_group = eventTable{prb}.with_stim == (g == 2);
        n_events(g) = sum(in_group);
        med_dur_ms(g) = median(eventTable{prb}.dur_ms(in_group));
        med_peak_rms(g) = median(eventTable{prb}.peak_rms(in_group));
        med_pr(g) = median(eventTable{prb}.pr(in_group));
    end
    stimSortSummary{prb} = table(group_names, n_events, med_dur_ms, med_peak_rms, med_pr, ...
        'VariableNames', {'group','n_events','med_dur_ms','med_peak_rms','med_pr'});
end


%% Data PreProcessing: Participation ratio in pre and post stim windows
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


%% Data Visualization: Events by stim group and pre vs post participation ratio
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
