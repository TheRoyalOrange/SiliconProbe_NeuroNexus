function oscillation_plotEvents(ProbeInfo, oscillations, probeSnippets)
% oscillation_plotEvents.m
%
% Description: Interactive viewer for probe-level events (prbEvent) - pages
%   through one prbEvent at a time, plotting every analyzed channel's
%   filtered signal for that event's window, in the same stacked/offset
%   shape as oshkosh_detect_funct.m's own "plot all events across channels"
%   figure (each channel plotted as trace(i) - offset*(i-1)), but showing
%   only the one event currently selected instead of a whole file's worth
%   of events at once. Left/right arrow keys move to the previous/next
%   event; up/down arrow keys cycle the color scheme used to color
%   channels, through two kinds of scheme:
%   (1) Shank (discrete/qualitative colormap: lines) -
%   every channel found in the same column of the same probe's
%   ProbeInfo.ProbeMaps entry gets the same color (columns are shanks per
%   ProbeInfo's own convention - see OpenEphys_BaseAnalysis.m).
%   (2) Depth (continuous colormap: copper) - within each shank, channels
%   are colored by their position along the shank (physical row order in
%   ProbeInfo.ProbeMaps, after excluding that probe's own
%   ProbeInfo.Ch_Remove channels and renumbering the remaining ones
%   starting at 1 - e.g. if a shank's first two physical rows are removed
%   channels, its next channel down becomes position 1), using ONE
%   colormap range shared across every shank (not normalized per shank) -
%   so color always progresses in strict sequence through the
%   depth-grouped plot order below, regardless of how many channels
%   survive removal on each shank. (A per-shank-normalized range was tried
%   first but let a shank with fewer surviving channels hit "full
%   brightness" at a shallower tier than one with more, breaking the
%   sequence when interleaved by depth across shanks.) Both shank
%   assignment and depth position are
%   looked up per raw channel ID (from probeSnippets.chanals), not by the
%   row position within the snippet data - those are not guaranteed to be
%   in ProbeMaps/shank order. The stacked plot's vertical (line) order
%   also depends on the active scheme: under Shank, lines keep the
%   natural chanals order (already grouped by shank); under Depth, lines
%   are reordered so channels are grouped by depth position first and
%   shank second - e.g. every shank's position-1 channel plots together,
%   then every shank's position-2 channel, and so on - to make same-depth
%   activity across shanks easy to compare. Only the plotting order
%   changes; each line's color is still tied to its own channel identity.
%   The x-axis shows actual recording time (ms, since fs=1kHz) rather than
%   a local 1:window_length index - traces, shading, and stim xlines are
%   all plotted against this same absolute time axis. Also marks stimulus
%   TTL activity on the plotted window: a red xline (LineWidth 1.5) at
%   each stim time that falls inside the currently displayed window, and
%   a grey shaded band behind the traces covering every timepoint up to
%   3000 samples (3s at 1kHz) AFTER (not before) any stim time in that
%   file - even a stim time outside the displayed window still shades the
%   nearby edge of the window it's close to, if within reach.
%   NOTE (design calls made without re-confirming, since these are cheap
%   to change later): When probeSnippets covers more than one probe, all
%   probes' events are flattened into a single left/right-browsable
%   sequence in probe-then-index order (all of probe 1's events, then all
%   of probe 2's, etc.) - not interleaved chronologically across probes.
%
% Inputs:
%   ProbeInfo (struct) - this animal's ProbeInfo struct (see
%     oshkosh_detect_funct.m for field details). Used here for .ProbeMaps
%     (cell array, 1 x probenum; each cell a 2D matrix of raw channel IDs
%     shaped like the physical probe, columns = shanks, rows = physical
%     position along the shank) and .Ch_Remove (cell array, 1 x probenum;
%     per-probe raw channel IDs excluded for the depth color scheme's
%     position numbering).
%   oscillations (struct) - first output of oshkosh_detect_funct.m. Used
%     fields: .stim_times (cell, 1 x Nfiles - per-file absolute stim
%     sample indices, or [] if none), .prbEvent_file_index and
%     .prbEvent_ext_timestamps (cell, 1 x length(probes) - used to find
%     which file a displayed prbEvent came from and its window's absolute
%     start sample, so stim_times can be converted to this plot's local
%     x-axis).
%   probeSnippets (struct) - third output of oshkosh_detect_funct.m.
%     Used fields: .chanals, .probes, .event_data (cell, 1 x
%     length(probes); event_data{p}{k} is a [numel(chanals) x
%     window_length] filtered-signal matrix for probe p's k-th prbEvent).
%
% Outputs:
%   (none) - opens an interactive figure window. Scroll with left/right
%   (event) and up/down (color palette) arrow keys; close the figure to
%   stop.
%
% Dependencies: Expects oscillations and probeSnippets to have been
%   generated together by the same oshkosh_detect_funct.m call (so
%   prbEvent indices/ordering line up between them), using this SAME
%   ProbeInfo (so probeSnippets.chanals are valid raw channel IDs to look
%   up in ProbeInfo.ProbeMaps).

if isempty(probeSnippets.probes)
    error('oscillation_plotEvents:NoProbes', 'probeSnippets has no probes.');
end

%% assign each channel (by raw channel ID, not row position) to a shank
chanals = probeSnippets.chanals;
nChan = numel(chanals);
shank_key = nan(nChan,2); % [probe index into ProbeMaps, column index within that probe's map]
for i = 1:nChan
    for pr = 1:numel(ProbeInfo.ProbeMaps)
        [~, col] = find(ProbeInfo.ProbeMaps{pr} == chanals(i));
        if ~isempty(col)
            shank_key(i,:) = [pr, col(1)];
            break
        end
    end
end
if any(isnan(shank_key(:)))
    error('oscillation_plotEvents:UnmappedChannel', ...
        'One or more channels in probeSnippets.chanals were not found in any ProbeInfo.ProbeMaps entry.');
end
[~, ~, shank_group] = unique(shank_key, 'rows'); % shank_group(i) = shank index (1..nShanks) for channel i
nShanks = max(shank_group);

%% position of each channel along its shank (depth), excluding that probe's Ch_Remove channels
depth_position = nan(nChan,1);
for i = 1:nChan
    pr = shank_key(i,1);
    col = shank_key(i,2);
    shank_col_ids = ProbeInfo.ProbeMaps{pr}(:,col); % physical row order = physical position along the shank
    kept_ids = shank_col_ids(~ismember(shank_col_ids, ProbeInfo.Ch_Remove{pr}));
    depth_position(i) = find(kept_ids == chanals(i), 1);
end

%% flatten every probe's events into one browsable sequence (probe-then-index order)
flat_pk = [];
for p = 1:numel(probeSnippets.probes)
    nEvents = numel(probeSnippets.event_data{p});
    flat_pk = [flat_pk; [repmat(p,nEvents,1), (1:nEvents)']];
end
if isempty(flat_pk)
    error('oscillation_plotEvents:NoEvents', 'probeSnippets has no detected probe events to plot.');
end

%% set up interactive figure
ud.oscillations = oscillations;
ud.probeSnippets = probeSnippets;
ud.shank_group = shank_group;
ud.nShanks = nShanks;
ud.depth_position = depth_position;
ud.flat_pk = flat_pk;
ud.nTotal = size(flat_pk,1);
ud.curIdx = 1;
ud.schemeMode = {'shank','depth'}; %parallel to schemeCmap
ud.schemeCmap = {'lines','copper'};
ud.cmapIdx = 1;

fig = figure('KeyPressFcn', @keypress, 'Name', 'oscillation_plotEvents');
ud.ax = axes('Parent', fig);
setappdata(fig, 'ud', ud);
redraw(fig);

end

function keypress(fig, evt)
    ud = getappdata(fig, 'ud');
    switch evt.Key
        case 'rightarrow'
            ud.curIdx = min(ud.curIdx + 1, ud.nTotal);
        case 'leftarrow'
            ud.curIdx = max(ud.curIdx - 1, 1);
        case 'uparrow'
            ud.cmapIdx = mod(ud.cmapIdx, numel(ud.schemeCmap)) + 1;
        case 'downarrow'
            ud.cmapIdx = mod(ud.cmapIdx - 2, numel(ud.schemeCmap)) + 1;
        otherwise
            return
    end
    setappdata(fig, 'ud', ud);
    redraw(fig);
end

function channel_colors = getChannelColors(ud)
    cmap_name = ud.schemeCmap{ud.cmapIdx};
    scheme_mode = ud.schemeMode{ud.cmapIdx};
    channel_colors = zeros(numel(ud.shank_group), 3);
    if strcmp(scheme_mode, 'shank')
        colors = feval(cmap_name, ud.nShanks);
        for i = 1:numel(ud.shank_group)
            channel_colors(i,:) = colors(ud.shank_group(i),:);
        end
    else %depth - one global colormap range across all shanks, so color always progresses
         %in sequence in the depth-grouped plot order regardless of how many
         %channels survive removal on each shank (a per-shank-normalized range
         %would let a shank with fewer channels hit "full brightness" at a
         %shallower tier than one with more, breaking the sequence)
        max_pos = max(ud.depth_position);
        colors = feval(cmap_name, max_pos);
        for i = 1:numel(ud.shank_group)
            channel_colors(i,:) = colors(ud.depth_position(i),:);
        end
    end
end

function plot_order = getPlotOrder(ud)
    scheme_mode = ud.schemeMode{ud.cmapIdx};
    nCh = numel(ud.shank_group);
    if strcmp(scheme_mode, 'shank')
        plot_order = (1:nCh)'; %natural chanals order - already grouped by shank
    else %depth - grouped by depth position first, shank second, so same-depth channels plot together across shanks
        [~, plot_order] = sortrows([ud.depth_position, ud.shank_group, (1:nCh)']);
    end
end

function [stim_abs, near_stim] = getStimOverlay(ud, abs_t)
    % stim_abs: stim times (absolute recording time, same units as abs_t)
    %   that fall inside this window - for the red xlines.
    % near_stim: logical, same size as abs_t, true wherever that
    %   timepoint is within 3000 samples AFTER (not before) ANY stim time
    %   in this file, whether or not that stim time itself falls inside
    %   the window - for shading.
    window_len = numel(abs_t);
    stim_abs = [];
    near_stim = false(1, window_len);

    if ~isfield(ud.oscillations, 'stim_times')
        return
    end

    p = ud.flat_pk(ud.curIdx,1);
    k = ud.flat_pk(ud.curIdx,2);
    file_idx = ud.oscillations.prbEvent_file_index{p}(k);
    file_stim_times = ud.oscillations.stim_times{file_idx};
    if isempty(file_stim_times)
        return
    end

    stim_abs = file_stim_times(file_stim_times >= abs_t(1) & file_stim_times <= abs_t(end));

    margin = 3000;
    % post-stim only: a stim time before the window can still shade into it
    % (its forward arm may reach [st, st+margin]), so filter on that reach.
    relevant = file_stim_times(file_stim_times >= abs_t(1)-margin & file_stim_times <= abs_t(end));
    for st = relevant(:)'
        lo = max(abs_t(1), st) - abs_t(1) + 1;
        hi = min(abs_t(end), st + margin) - abs_t(1) + 1;
        near_stim(lo:hi) = true;
    end
end

function redraw(fig)
    ud = getappdata(fig, 'ud');
    p = ud.flat_pk(ud.curIdx,1);
    k = ud.flat_pk(ud.curIdx,2);
    snippet = ud.probeSnippets.event_data{p}{k};
    nCh = size(snippet,1);
    cmap_name = ud.schemeCmap{ud.cmapIdx};
    scheme_mode = ud.schemeMode{ud.cmapIdx};
    channel_colors = getChannelColors(ud);
    plot_order = getPlotOrder(ud);
    offset_step = 500;
    window_len = size(snippet,2);
    file_idx = ud.oscillations.prbEvent_file_index{p}(k);
    ext_start = ud.oscillations.prbEvent_ext_timestamps{p}(k,1);
    abs_t = ext_start : (ext_start + window_len - 1); % actual recording timestamps (ms, since fs=1kHz)
    [stim_abs, near_stim] = getStimOverlay(ud, abs_t);

    cla(ud.ax);
    hold(ud.ax, 'on');

    %shading behind the traces: grey band near any stim time (in or out of this window)
    if any(near_stim)
        d = diff([0, near_stim, 0]);
        run_starts = find(d == 1);
        run_ends = find(d == -1) - 1;
        yhi = max(snippet(:));
        ylo = min(snippet(:)) - offset_step*(nCh-1);
        pad = 0.05*(yhi - ylo + eps);
        for r = 1:numel(run_starts)
            patch(ud.ax, abs_t([run_starts(r) run_ends(r) run_ends(r) run_starts(r)]), ...
                [ylo-pad ylo-pad yhi+pad yhi+pad], [0.7 0.7 0.7], ...
                'FaceAlpha', 0.4, 'EdgeColor', 'none');
        end
    end

    for slot = 1:nCh
        i = plot_order(slot);
        plot(ud.ax, abs_t, snippet(i,:) - offset_step*(slot-1), 'Color', channel_colors(i,:));
    end

    %red xline for each stim time that actually falls inside this window
    for sa = stim_abs(:)'
        xline(ud.ax, sa, 'Color', 'r', 'LineWidth', 1.5);
    end

    hold(ud.ax, 'off');
    xlabel(ud.ax, 'Time in Recording (ms)');
    title(ud.ax, sprintf('Probe %d - Recording %d: Event %d/%d  (overall %d/%d)  -  colors: %s (%s)', ...
        ud.probeSnippets.probes(p), file_idx, k, numel(ud.probeSnippets.event_data{p}), ud.curIdx, ud.nTotal, scheme_mode, cmap_name), ...
        'Interpreter', 'none');
end
