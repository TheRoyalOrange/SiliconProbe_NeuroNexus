% ConditionalTrialRemove_callable.m
%
% Description: Callable, loop-friendly trial-curation utility. For one
%   animal/condition/channel, checks whether a NAMED trial-exclusion
%   condition (tr_conditional) already has a recorded trial set in this
%   file's own tr_remove_conditional table (table, columns Name/
%   TrialIdx, stored inside the LFP/spiking/TF/CSD result files
%   themselves - see Outputs). If a match already exists, does nothing
%   and returns immediately (no figures, no dialogs, no save) - this is
%   what makes it fast to call in a loop over many files, most of which
%   already have this condition recorded. If no match exists, tells the
%   user via a dialog, then lets them visually inspect LFP and
%   spike-raster traces for one channel of an already-preprocessed
%   animal/condition (output of OpenEphys_BaseAnalysis.m or a variant),
%   pick trials to exclude, and records that choice under the
%   tr_conditional name. Does NOT touch tr_keep, and does NOT modify the
%   LFP/spiking/TF/CSD result files beyond the tr_remove_conditional
%   write described in Outputs and one other exception: if tr_remove is
%   missing or empty in spiking_results.mat (old-format data that
%   predates trial-removal tracking), it is migrated in place -
%   defaulted to all-zeros (no trials removed) and written back, via
%   save(...,'-append'), to the LFP/spiking/TF/CSD result files for this
%   animal/condition, so all four stay in sync. The recorded exclusion
%   set itself is applied later, at analysis time, by whichever script
%   chooses to combine it with tr_remove into a tr_remove_local mask.
%   Separately, every call (whether it skips via a match or records a
%   new selection) also checks/updates a single master tracking file,
%   <save_directory>\tr_conditional_master.mat (table, columns Name/
%   Animal/Filename - no TrialIdx, that stays in the per-file tables) -
%   an index of which tr_conditional names exist and which animal/
%   filename pairs each has been applied to, so that can be scanned on
%   its own without opening every per-recording result file. See
%   updateTrConditionalMaster (local function, bottom of this file).
%
% Inputs:
%   animal (string) - animal name matching the folder under
%     save_directory that holds this recording's result files.
%   condition (string) - condition name only (e.g. 'whisker'), NOT
%     animal-prefixed. Internally combined with animal as
%     full_filename = [animal '-' condition] (e.g.
%     '20260423-p12-whisker'), matching the <animal>-<condition> prefix
%     that OpenEphys_BaseAnalysis.m (or a variant) actually names the
%     LFP/spiking/figure files with.
%   tr_conditional (string) - REQUIRED. Name of the trial-exclusion
%     condition to look up/record for this file (e.g.
%     'whisker_twitch_artifact'). Not optional and has no global/
%     persistent fallback - the caller is expected to already know what
%     they want to select before starting (e.g. set once in their own
%     loop/script and pass it in on every call), since MATLAB doesn't
%     allow a global or persistent variable to share a name with a
%     function parameter, which ruled out a "prompt once, remember
%     across calls" fallback under this same name.
%   fs (int, optional) - sampling rate of the original recording, in Hz
%     (used to bin spikes into 1 ms bins). Defaults to 30000.
%   ProbeInfo (struct, optional) - This animal's ProbeInfo struct (as
%     loaded from <animal>-ProbeInfo.mat). If omitted or left empty
%     ([]), loaded automatically from
%     <save_directory>\ProbeInfo\<animal>-ProbeInfo.mat. Used only for
%     .Areas (to locate the LFP overview figures below) - no longer
%     carries tr_remove_conditional (moved to per-file storage; see
%     Outputs).
%   save_directory (string, optional) - base path containing the LFP/
%     Spiking/ProbeInfo subfolders. Defaults to
%     'E:\Roy\Processed Silicon Probe Data'.
%   master_filename (string, optional) - filename (not a full path) of
%     the master tracking file within save_directory - see
%     updateTrConditionalMaster below. Defaults to
%     'tr_conditional_master.mat'. Exists so a caller can point at a
%     throwaway/test master file (e.g. when dry-running against test
%     data) without touching the real one - existing callers that omit
%     it are unaffected.
%   tr_keep (int vector, trials), tr_remove (0/1 vector, same length as
%     tr_keep), tr_remove_conditional (table, columns Name/TrialIdx;
%     defaults to an empty table if missing from the file. TrialIdx may
%     be NaN for a given Name - that means this condition was evaluated
%     for this file and the user confirmed zero trials removed (keep
%     all), not a literal trial index; any code that later reads this
%     table to build a tr_remove_local mask must treat a NaN TrialIdx as
%     a no-op rather than an index) - loaded automatically, up front, in
%     one lightweight call, from
%     <save_directory>\Spiking\<animal>\<full_filename>-spiking_results.mat.
%     If tr_remove is missing or empty (old-format file), it is
%     migrated to zeros(1,length(tr_keep)) and written back - see
%     Description. This initial load deliberately excludes
%     stim_spike_stimchunks (the heavy 30kHz array) - see chan below for
%     when/how that's read.
%   stim_lfp_stimchunks - loaded automatically from
%     <save_directory>\LFP\<animal>\<full_filename>-LFP.mat. 3D array of
%     lfp data split into trials, downsampled to 1kHz. size[trials x
%     trial length x channels]. Loaded in full (unlike
%     stim_spike_stimchunks below) because this file is not saved with
%     '-v7.3', so a matfile partial-channel read would still have to
%     pull the whole array off disk anyway - no benefit to deferring it.
%   do_it_live (0/1, optional) - defaults to 0. 1 = skip the activity
%     label pick (and its no-labels dialog) and go straight to the old
%     channel pick (LFP overview figs + typed channel number) - see chan
%     below. ConditionalTrialRemove.m passes 1 for animals whose ProbeInfo
%     has no labels, after the user already said "do it live" in its
%     pre-run check, so they aren't asked again for every file.
%   ProbeInfo activity labels (optional, used to pick chan) - read from
%     the ProbeInfo struct above. Made by ChannelShankLabelling_StimData.m
%     / ChannelShankLabelling_SponData.m:
%       .ChanLabels (cell, 1 x probenum) - each a struct, one field per
%         label name; field = double [rows x shanks] channel scores in the
%         ProbeMaps layout, NaN for Ch_Remove channels
%       .ShankLabels (cell, 1 x probenum) - same, field = double
%         [1 x shanks] shank scores (optional - the line plot is just left
%         empty if a label has no shank scores)
%       .ProbeMaps (cell, 1 x probenum) - [rows x shanks] raw channel IDs
%       .Ch_Remove (cell, 1 x probenum) - raw channel IDs to exclude
%   chan (int) - NOT a function input. Chosen after the LFP data is
%     loaded, so the caller doesn't need to know the channel count in
%     advance:
%     - If ProbeInfo.ChanLabels has any labels, the user picks one from a
%       list, then gets a channel map figure (pickChanFromLabelMap, local
%       function): shank score line plot on top, channel score heatmap in
%       the ProbeMaps layout below, Ch_Remove channels black, channel IDs
%       overlaid. The user clicks the channel they want, then confirms.
%     - If there are no labels (or the user cancels the label list, or
%       closes the map without picking), a two-button dialog
%       (twoButtonDialog) says so. "Fuck it, we'll do it live" falls
%       back to the original method: open the saved LFP overview
%       figure(s) for this animal/condition (one per probe/region, from
%       <save_directory>\Figures_BaseAnalysis\<animal>\
%       <full_filename>-<ProbeInfo.Areas{prb}>_Probe-lfp_results.fig)
%       and type the channel number into an inputdlg (bounded 1 to the
%       number of channels in stim_lfp_stimchunks). "Cancel. I'll be
%       back" returns right away with runCancelled = true, saving nothing.
%     INFERRED DATA CONTRACT (new with the channel map): a clicked map
%     cell gives chan = ProbeInfo.ProbeMaps{prb}(r,c), the raw channel
%     ID. This assumes the 3rd dim of stim_lfp_stimchunks /
%     stim_spike_stimchunks is indexed by raw channel ID, the same way
%     ChannelShankLabelling_StimData.m indexes it to make the labels.
%     Once chosen,
%     stim_spike_stimchunks is read for this channel only, via a
%     matfile object indexed as stim_spike_stimchunks(:,:,chan), so only
%     one channel's data is ever held in memory afterward (checked via
%     checkcode + a synthetic-data correctness test - bit-identical to
%     the old all-channel-then-slice approach). NOT a disk I/O win,
%     though: checked with h5info against a real saved file, MATLAB's
%     default '-v7.3' auto-chunking on this array picks a chunk shape
%     spanning ALL trials and ALL channels (e.g. [40 6 32] for a
%     [40 trials x 105000 samples x 32 channels] array - just a thin
%     slice of samples per chunk), so every chunk has to be touched
%     regardless of which channel is requested; a timed comparison
%     against a full load showed no measurable speed difference. Fixing
%     that would require re-chunking OpenEphys_BaseAnalysis.m's save
%     calls and re-saving existing files - out of scope here. Binning
%     below operates on just this channel's slice. INFERRED DATA
%     CONTRACT (pre-existing, not introduced by this change - the old
%     all-channel binning loop and this per-channel read both assume
%     it): chan
%     indexes the same channel in both stim_lfp_stimchunks and
%     stim_spike_stimchunks, i.e. the two arrays have matching channel
%     count and ordering.
%
% Outputs:
%   ProbeInfo (struct) - the input struct, returned largely unchanged
%     (still useful to the caller for .Areas etc. without reloading it
%     themselves in a loop). No longer carries tr_remove_conditional.
%   runCancelled (logical scalar) - true only if the user picked
%     "Cancel. I'll be back" at the channel pick (see chan above), or
%     "Abort. I'll be back" at the missing-result-file check (see
%     fileSkipped below); nothing is saved for this file in that case and
%     the caller should stop its loop. false otherwise, including the
%     "match found, skip" path.
%   fileSkipped (logical scalar) - existence of this animal/condition's
%     four result files (spiking_results, LFP, CSD_results, _TF_results -
%     see tr_remove_conditional below) is checked first thing, before any
%     load or dialog. If the spiking or LFP file is missing (the trial
%     plots need both), a non-modal dialog lists the missing file(s) and
%     offers "Skip this file, keep going" (returns right away, fileSkipped
%     = true, nothing saved or added to the master index) or "Abort. I'll
%     be back" (returns right away, runCancelled = true; closing the box
%     counts as abort). false otherwise.
%   missingsaves (cell array of strings, 1 x missing files) - if only the
%     CSD and/or TF file is missing, the dialog instead offers "Skip saving
%     to those, keep going": everything else still runs (channel pick,
%     trial pick, master index update) and the tr_remove migration /
%     tr_remove_conditional saves just skip the missing file(s), whose
%     names (no folder) come back here. Abort works as above. Empty if all
%     four exist.
%   (side effect, not a return value) tr_remove_conditional - if no
%     match was found for tr_conditional, the newly chosen trial set is
%     appended (as Name/TrialIdx rows) and written, identically, into
%     all four of this animal/condition's result files:
%     <full_filename>-spiking_results.mat, <full_filename>-LFP.mat,
%     <full_filename>-CSD_results.mat, and
%     <full_filename>_TF_results.mat (note: no dash before
%     "_TF_results", unlike the other three - matches how
%     OpenEphys_BaseAnalysis.m actually names that file). Also logged to
%     the console.
%   (side effect, not a return value) tr_conditional_master - the
%     Name/Animal/Filename row for this call is added if not already
%     present, in <save_directory>\<master_filename> (see Inputs -
%     defaults to tr_conditional_master.mat). Updated on every call,
%     including the "match found, skip" path. Also logged to the
%     console.
%
% Dependencies: Expects OpenEphys_BaseAnalysis.m (or the _Bundled /
%   _MixedTrials variant) to have already been run for this
%   animal/condition, producing the LFP.mat and spiking_results.mat
%   files this function reads, and <animal>-ProbeInfo.mat to already
%   exist. Also expects ProbeInfo.Areas (cell array of probe/region
%   names) to be populated, and (non-fatal if missing - a warning is
%   issued instead) the per-probe lfp_results.fig files under
%   Figures_BaseAnalysis to have been saved for this animal/condition (only
%   opened on the "do it live" fallback). Optionally,
%   ChannelShankLabelling_StimData.m / ChannelShankLabelling_SponData.m
%   run first so ProbeInfo has activity labels for the channel map.

function [ProbeInfo, runCancelled, fileSkipped, missingsaves] = ConditionalTrialRemove_callable(animal, condition, tr_conditional, fs, ProbeInfo, save_directory, master_filename, do_it_live)

runCancelled = false; %flips to true only if the user hits "Cancel. I'll be back" at the channel pick, or "Abort" at the missing-file check
fileSkipped = false; %flips to true only if the spiking or LFP file is missing and the user says skip this file
missingsaves = {}; %names of missing result files that didn't get saved into (filled in by the check below)

if nargin < 8 || isempty(do_it_live)
    do_it_live = 0; %default: try the activity label channel map first
end

if nargin < 7 || isempty(master_filename)
    master_filename = 'tr_conditional_master.mat';
end

if nargin < 6 || isempty(save_directory)
    save_directory = 'E:\Roy\Processed Silicon Probe Data';
end

if nargin < 5 || isempty(ProbeInfo)
    probeinfo_dat = load(fullfile([save_directory '\ProbeInfo\' animal '-ProbeInfo.mat']), 'ProbeInfo');
    ProbeInfo = probeinfo_dat.ProbeInfo;
end

if ~isstruct(ProbeInfo) || ~isscalar(ProbeInfo)
    error('ConditionalTrialRemove_callable:InvalidProbeInfo', ...
        ['ProbeInfo must be a scalar struct, got a %s. This usually means a stale call ' ...
         'still passes chan as the 3rd argument - the current signature is ' ...
         'ConditionalTrialRemove_callable(animal, condition, tr_conditional, fs, ProbeInfo, save_directory), ' ...
         'with chan now chosen via dialog instead of as an input.'], class(ProbeInfo));
end

if nargin < 4 || isempty(fs)
    fs = 30000;
end

if ~isnumeric(fs) || ~isscalar(fs)
    error('ConditionalTrialRemove_callable:InvalidFs', ...
        ['fs must be a numeric scalar, got a %s. This usually means fs was omitted from a call ' ...
         'that still supplied later arguments, letting one of them (often ProbeInfo) land in the ' ...
         'fs slot instead - the current signature is ConditionalTrialRemove_callable(animal, ' ...
         'condition, tr_conditional, fs, ProbeInfo, save_directory). Pass [] for fs to use its ' ...
         'default (30000) while still supplying later arguments.'], class(fs));
end

full_filename = [animal '-' condition]; % matches the <animal>-<condition> prefix OpenEphys_BaseAnalysis.m actually writes files under
spikingFile = fullfile([save_directory '\Spiking\' animal '\' full_filename '-spiking_results.mat']);

%% make sure all four result files exist before anything else (the pick gets saved into all four, so no point picking otherwise)
resultFiles = {spikingFile, ...
    fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']), ...
    fullfile([save_directory '\CSD\' animal '\' full_filename '-CSD_results.mat']), ...
    fullfile([save_directory '\TF\' animal '\' full_filename '_TF_results.mat'])}; %cell {1 x 4} spiking, LFP, CSD, TF
fileexists = false(1,numel(resultFiles)); %[1 x 4] 1 = file is there, gets saved into; 0 = missing, save skipped
%missingsaves (made empty at the top) gets the file names (no folder) of whichever ones aren't there
for i = 1:numel(resultFiles)
    fileexists(i) = isfile(resultFiles{i});
    if fileexists(i) == false
        [~, nm, ext] = fileparts(resultFiles{i});
        missingsaves{end+1} = [nm ext]; %#ok<AGROW> at most 4
    end
end

%spiking or LFP missing? can't plot trials without them - skip this whole file, or abort the run
if fileexists(1) == false || fileexists(2) == false
    missingMsg = sprintf(['Can''t pick trials for %s - no file found for:\n\n%s\n\n' ...
        '(need the spiking + LFP files to plot.) Skip this file and keep going, or abort the whole run?'], full_filename, strjoin(missingsaves, '\n'));
    missingAnswer = twoButtonDialog(missingMsg, 'Result file(s) missing', 'Skip this file, keep going', 'Abort. I''ll be back');
    if strcmp(missingAnswer, 'Skip this file, keep going')
        fileSkipped = true; %caller adds it to its skipped list for the post-run report
    else %abort button, or closed the box
        disp(['ConditionalTrialRemove_callable: run aborted at ' full_filename ' (result file(s) missing) - nothing saved for this file.']);
        runCancelled = true;
    end
    return
end
%(only CSD and/or TF missing is asked about further down, after the "match found" check, so already-done files don't ask)

%% light load: just enough to decide whether to skip and to migrate tr_remove if needed - deliberately excludes
%% stim_spike_stimchunks (the heavy 30kHz array), so a "match found, skip" call never touches it. See below for
%% where the chosen channel's spike data is read on its own, only once we know we're not skipping.
spiking_dat = load(spikingFile, 'tr_keep','tr_remove','tr_remove_conditional');
tr_keep = spiking_dat.tr_keep;

if ~isfield(spiking_dat,'tr_remove_conditional')
    tr_remove_conditional = table('Size',[0,2], 'VariableTypes',{'string','double'}, 'VariableNames',{'Name','TrialIdx'});
else
    tr_remove_conditional = spiking_dat.tr_remove_conditional;
end

%% if this named condition already has a recorded trial set for this file, skip everything below (loop-friendly: no dialogs, no figures, no save)
if any(strcmp(tr_remove_conditional.Name, tr_conditional))
    fprintf('Match found for ''%s'' in %s - skipping.\n', tr_conditional, full_filename);
    updateTrConditionalMaster(save_directory, master_filename, tr_conditional, animal, full_filename);
    missingsaves = {}; %nothing was going to be saved anyway, so don't report these as skipped
    return
end

%only CSD and/or TF missing? still pick trials, just don't save into those - or abort the run
if ~isempty(missingsaves)
    missingMsg = sprintf(['No file found for:\n\n%s\n\n''%s'' won''t be saved into these for %s. ' ...
        'Skip saving to those and keep going, or abort the whole run?'], strjoin(missingsaves, '\n'), tr_conditional, full_filename);
    missingAnswer = twoButtonDialog(missingMsg, 'Result file(s) missing', 'Skip saving to those, keep going', 'Abort. I''ll be back');
    if ~strcmp(missingAnswer, 'Skip saving to those, keep going') %abort button, or closed the box
        disp(['ConditionalTrialRemove_callable: run aborted at ' full_filename ' (result file(s) missing) - nothing saved for this file.']);
        missingsaves = {};
        runCancelled = true;
        return
    end
end
% uiwait blocks script execution until the user clicks OK (or closes the box);
% msgbox itself stays non-modal (WindowStyle 'normal'), so other figures are still clickable while it waits.
uiwait(msgbox(sprintf('No match found for ''%s'' in %s. Please select trials to remove.', tr_conditional, full_filename)));

%% migrate tr_remove from the old format (missing/empty field) to the new format (0/1 vector matching tr_keep), writing the fix back to every related result file
if ~isfield(spiking_dat,'tr_remove') || isempty(spiking_dat.tr_remove)
    tr_remove = zeros(1,length(tr_keep));
    % save(spikingFile,'tr_remove','-append');
    % save(fullfile([save_directory '\CSD\' animal '\' full_filename '-CSD_results.mat']),'tr_remove','-append');
    % save(fullfile([save_directory '\TF\' animal '\' full_filename '_TF_results.mat']),'tr_remove','-append');
    % save(fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']),'tr_remove','-append');
    %^replaced by the loop below: only saves into the result files that exist (missing CSD/TF get skipped, see check up top)
    for i = find(fileexists)
        save(resultFiles{i},'tr_remove','-append');
    end
else
    tr_remove = spiking_dat.tr_remove;
end
lfp_dat = load(fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']), 'stim_lfp_stimchunks');
stim_lfp_stimchunks = lfp_dat.stim_lfp_stimchunks;

tr_keep_local = tr_keep(~logical(tr_remove)); %only offer trials not already permanently excluded

%% pick which channel to inspect, off this animal's probe activity labels (ProbeInfo.ChanLabels) if it has any
% Every dialog in this function is kept non-modal where MATLAB allows it
% (WindowStyle 'normal' on inputdlg calls; questdlg has no non-modal
% option at all, so its Yes/No use below is replaced with the
% nonModalConfirm() local function), so the user can always click into
% the LFP overview figure(s) or trial plots while deciding. listdlg is
% modal, but it's only used for the label pick, before any figure opens.

%what activity labels does this animal have? (every label name found on any probe)
labelnames = {};
if isfield(ProbeInfo,'ChanLabels')
    for prb = 1:numel(ProbeInfo.ChanLabels)
        if isstruct(ProbeInfo.ChanLabels{prb})
            labelnames = union(labelnames, fieldnames(ProbeInfo.ChanLabels{prb}), 'stable'); %cell {label names x 1}
        end
    end
end

chan = []; %stays empty until a channel is picked off the label map
if do_it_live == 0 %do_it_live = 1 means the caller already chose the old way - skip the label pick entirely
    if isempty(labelnames)
        fallbackMsg = sprintf(['No probe activity labels exist for %s.\n\nRecommend running ChannelShankLabelling_StimData.m ' ...
            'and/or ChannelShankLabelling_SponData.m first to make them.'], animal);
    else
        %which label should color the channel map?
        [labelIdx, labelOK] = listdlg('ListString', labelnames, 'SelectionMode', 'single', 'ListSize', [300 150], ...
            'Name', 'Pick an activity label', 'PromptString', sprintf('Which activity label for %s?', full_filename));
        if labelOK == 0
            fallbackMsg = sprintf('No activity label chosen for %s.\n\nPick a channel the old way, or bail on the whole run?', full_filename);
        else
            %plot the label as a probe-shaped channel map and let the user click their channel
            chan = pickChanFromLabelMap(ProbeInfo, labelnames{labelIdx}, full_filename);
            fallbackMsg = sprintf('Channel map closed without picking a channel for %s.\n\nPick a channel the old way, or bail on the whole run?', full_filename);
        end
    end
end

%no channel from the label map? do it live (the original way) or cancel the whole run
if isempty(chan)
    if do_it_live == 0 %already chosen up front if do_it_live = 1, so don't ask again
        fallbackAnswer = twoButtonDialog(fallbackMsg, 'No channel picked yet', 'Fuck it, we''ll do it live', 'Cancel. I''ll be back');
        if strcmp(fallbackAnswer, 'Cancel. I''ll be back')
            disp(['ConditionalTrialRemove_callable: run cancelled at ' full_filename ' - nothing saved for this file.']);
            runCancelled = true;
            return
        end
    end

    %open the saved LFP overview figure(s) for this animal/condition (one per probe/region), for channel reference
    for prb = 1:numel(ProbeInfo.Areas)
        %lfpFigFile = fullfile([save_directory '\AnimalFigures\' animal '\' full_filename '-' ProbeInfo.Areas{prb} '_Probe-lfp_results.fig']); %old folder name, BaseAnalysis now saves to Figures_BaseAnalysis
        lfpFigFile = fullfile([save_directory '\Figures_BaseAnalysis\' animal '\' full_filename '-' ProbeInfo.Areas{prb} '_Probe-lfp_results.fig']);
        if isfile(lfpFigFile)
            openfig(lfpFigFile);
        else
            warning('ConditionalTrialRemove_callable:LFPFigNotFound', 'LFP overview figure not found for %s, probe ''%s'': %s', full_filename, ProbeInfo.Areas{prb}, lfpFigFile);
        end
    end

    %ask which channel to inspect
    nchan = size(stim_lfp_stimchunks,3);
    chanAnswer = inputdlg({sprintf('Which channel would you like to inspect (1-%d)?', nchan)}, 'Select channel', [1 50], {''}, struct('WindowStyle','normal'));
    chan = str2double(chanAnswer{1});
end

%% read only the chosen channel's spike data (keeps just 1 channel in memory afterward instead of all
%% of them; NOT a disk I/O win here - see header doc re: this array's HDF5 chunk layout)
spikingMat = matfile(spikingFile);
chspikes = squeeze(spikingMat.stim_spike_stimchunks(:,:,chan));

%% bin that channel's spikes into 1ms bins for plotting
spikes_ms = [];
for batch = 1:size(chspikes,2)/(fs/1000)
    spikebatchi = sum(chspikes(:,1+((fs/1000)*(batch-1)):(fs/1000)+((fs/1000)*(batch-1))),2)>0;
    spikes_ms(:,batch) = spikebatchi;
end

%% pick trials, with a confirm/redo loop
time = [-500:3000];
confirmed = false;
while ~confirmed
    figure();
    for tr = 1:length(tr_keep_local)
        subplot(1,2,1)
        hold on
        colormap(sky(8))
        plot(time,squeeze(stim_lfp_stimchunks(tr_keep_local(tr),4500:8000,chan))'-(400*(tr-1)))
        xline(0)
        subplot(1,2,2)
        hold on
        plot(time,spikes_ms(tr_keep_local(tr),4500:8000)'-(2*(tr-1)))
        xline(0)
    end

    promt = {sprintf('Which of these trials (plot row number, 1-%d) would you like to remove for condition ''%s''?', length(tr_keep_local), tr_conditional)};
    removed_rownum = inputdlg(promt, 'Conditional trial removal', [1 50], {''}, struct('WindowStyle','normal'));
    removed_rownum = str2num(removed_rownum{1});
    chosen_trialidx = tr_keep_local(removed_rownum);

    quest = sprintf('Remove %d trial(s) under condition ''%s''?', numel(chosen_trialidx), tr_conditional);
    answer = nonModalConfirm(quest, 'Confirm conditional removal');
    close all
    if strcmp(answer, 'Yes')
        confirmed = true;
    end
end

%% record the new selection (existingRows was already confirmed empty by the match check above)
% zero trials chosen (keep all) still gets one row, TrialIdx = NaN, so this
% condition/file pair is recorded as evaluated and isn't re-prompted next time
if isempty(chosen_trialidx)
    newRows = table(string(tr_conditional), NaN, 'VariableNames', {'Name','TrialIdx'});
else
    newRows = table(repmat(string(tr_conditional),numel(chosen_trialidx),1), chosen_trialidx(:), ...
        'VariableNames', {'Name','TrialIdx'});
end
tr_remove_conditional = [tr_remove_conditional; newRows];

%% save the new selection into all four result files, and log what was written
fprintf('\n--- ConditionalTrialRemove_callable: rows added for %s, condition ''%s'' ---\n', full_filename, tr_conditional);
disp(newRows);

% save(fullfile([save_directory '\Spiking\' animal '\' full_filename '-spiking_results.mat']), 'tr_remove_conditional', '-append');
% save(fullfile([save_directory '\LFP\' animal '\' full_filename '-LFP.mat']), 'tr_remove_conditional', '-append');
% save(fullfile([save_directory '\CSD\' animal '\' full_filename '-CSD_results.mat']), 'tr_remove_conditional', '-append');
% save(fullfile([save_directory '\TF\' animal '\' full_filename '_TF_results.mat']), 'tr_remove_conditional', '-append');
%^replaced by the loop below: only saves into the result files that exist (missing CSD/TF get skipped, see check up top)
for i = find(fileexists)
    save(resultFiles{i}, 'tr_remove_conditional', '-append');
end

updateTrConditionalMaster(save_directory, master_filename, tr_conditional, animal, full_filename);

end

function answer = nonModalConfirm(question, dlgtitle)
% Minimal Yes / No, redo selection confirm dialog. questdlg hard-codes
% WindowStyle to 'modal' with no public override, so it can't be made
% non-modal directly; this reimplements just enough of it (same two
% button labels/behavior) as an ordinary, non-modal figure so the user
% can still click into the trial-plot figure while deciding.
answer = 'No, redo selection';
fig = figure('Name', dlgtitle, 'NumberTitle', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'WindowStyle', 'normal', 'Resize', 'off', ...
    'Position', [400 400 340 110]);
uicontrol(fig, 'Style', 'text', 'String', question, 'Units', 'normalized', ...
    'Position', [0.05 0.5 0.9 0.4], 'HorizontalAlignment', 'center');
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Yes', 'Units', 'normalized', ...
    'Position', [0.1 0.15 0.3 0.3], 'Callback', @(~,~) respond('Yes'));
uicontrol(fig, 'Style', 'pushbutton', 'String', 'No, redo selection', 'Units', 'normalized', ...
    'Position', [0.55 0.15 0.35 0.3], 'Callback', @(~,~) respond('No, redo selection'));
uiwait(fig);
    function respond(val)
        answer = val;
        uiresume(fig);
        close(fig);
    end
end

function answer = twoButtonDialog(question, dlgtitle, btn1, btn2)
% Two-button, non-modal dialog (was chanPickerFallback, now shared with the
% missing-result-file check). Used when no channel came off the activity
% label map (btn1 = do it live, btn2 = cancel) and when a result file is
% missing (btn1 = skip it, btn2 = abort). Same build as nonModalConfirm
% above. Closing the box counts as btn2 (the cancel/abort one), so nothing
% gets saved by accident.
answer = btn2;
fig = figure('Name', dlgtitle, 'NumberTitle', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'WindowStyle', 'normal', 'Resize', 'off', ...
    'Position', [400 400 440 170]);
uicontrol(fig, 'Style', 'text', 'String', question, 'Units', 'normalized', ...
    'Position', [0.05 0.42 0.9 0.53], 'HorizontalAlignment', 'center');
uicontrol(fig, 'Style', 'pushbutton', 'String', btn1, 'Units', 'normalized', ...
    'Position', [0.05 0.1 0.42 0.25], 'Callback', @(~,~) respond(btn1));
uicontrol(fig, 'Style', 'pushbutton', 'String', btn2, 'Units', 'normalized', ...
    'Position', [0.53 0.1 0.42 0.25], 'Callback', @(~,~) respond(btn2));
uiwait(fig);
    function respond(val)
        answer = val;
        uiresume(fig);
        close(fig);
    end
end

function chan = pickChanFromLabelMap(ProbeInfo, labelname, full_filename)
% Plots one activity label as a probe-shaped channel map and lets the user
% click the channel to inspect. One column per probe that has this label:
% top = shank score line plot (ProbeInfo.ShankLabels), bottom = channel
% score heatmap in the ProbeMaps layout (ProbeInfo.ChanLabels), Ch_Remove
% channels black with an 'X', every other cell labelled with its channel
% ID. Same look as ChannelShankLabelling_StimData.m's plot. Each click is
% confirmed with nonModalConfirm; "No, redo selection" waits for another
% click. Returns the picked raw channel ID (pmap(r,c)), or [] if the user
% closes the figure without picking.
chan = [];

%which probes have this label?
prbs = [];
for prb = 1:numel(ProbeInfo.ChanLabels)
    if isstruct(ProbeInfo.ChanLabels{prb}) && isfield(ProbeInfo.ChanLabels{prb}, labelname)
        prbs(end+1) = prb; %#ok<AGROW> only a handful of probes
    end
end

seqmap = interp1([0 0.5 1],[0.95 0.98 0.93; 0.55 0.82 0.40; 0.18 0.60 0.22],linspace(0,1,256)); %same pale mint -> green colormap as ChannelShankLabelling_StimData.m
fig = figure('Name', [full_filename ' - pick a channel'], 'NumberTitle', 'off', 'Color', 'w');
tl = tiledlayout(fig, 4, numel(prbs), 'TileSpacing', 'compact');
title(tl, [full_filename ' - ' labelname ': click the channel to inspect'], 'Interpreter', 'none');

for p = 1:numel(prbs)
    prb = prbs(p);
    pmap = ProbeInfo.ProbeMaps{prb}; %[rows x shanks] raw channel IDs
    [nrow,nshnk] = size(pmap);
    removed = ismember(pmap, ProbeInfo.Ch_Remove{prb});
    chScoreMap = ProbeInfo.ChanLabels{prb}.(labelname); %[rows x shanks], ProbeMaps layout
    chScoreMap(removed) = NaN; %in case Ch_Remove changed after the label was made

    %top: shank scores (if this label has them)
    ax1 = nexttile(tl, p);
    if numel(ProbeInfo.ShankLabels) >= prb && isstruct(ProbeInfo.ShankLabels{prb}) && isfield(ProbeInfo.ShankLabels{prb}, labelname)
        plot(ax1, 1:nshnk, ProbeInfo.ShankLabels{prb}.(labelname), '-o', 'LineWidth', 2, 'MarkerSize', 8, 'Color', seqmap(end,:), 'MarkerFaceColor', seqmap(end,:));
    end
    xlim(ax1, [0.5 nshnk+0.5]);
    xticks(ax1, 1:nshnk); xticklabels(ax1, {});
    ylabel(ax1, {'Shank score','(mean of channels)'});
    title(ax1, ['Probe ' num2str(prb)]);
    grid(ax1, 'on'); ax1.GridAlpha = 0.15; box(ax1, 'off');

    %bottom: clickable channel heatmap (NaN cells are transparent over a black axes, so removed = black)
    ax2 = nexttile(tl, numel(prbs)+p, [3 1]);
    img = imagesc(ax2, 1:nshnk, 1:nrow, chScoreMap, 'AlphaData', ~isnan(chScoreMap));
    img.PickableParts = 'all'; %transparent (removed) cells still catch the click, so we can say "pick another"
    img.ButtonDownFcn = @(src,~) chanMapClick(src, pmap, chScoreMap, fig);
    ax2.Color = 'k';
    colormap(ax2, seqmap);
    cb = colorbar(ax2);
    if numel(prbs) == 1
        cb.Layout.Tile = 'east'; %keep colorbar outside both tiles so the shank columns stay aligned
    end
    cb.Label.String = ['Channel score (' labelname ')'];
    cb.Label.Interpreter = 'none';
    xticks(ax2, 1:nshnk); yticks(ax2, 1:nrow);
    xlabel(ax2, 'Shank (ProbeMaps column)'); ylabel(ax2, 'ProbeMaps row');
    linkaxes([ax1 ax2], 'x');

    %cell labels: channel ID, or 'X' for removed; black on green, white on black. Not clickable, so clicks reach the image underneath
    for r = 1:nrow
        for c = 1:nshnk
            if removed(r,c)
                txt = 'X';
            else
                txt = num2str(pmap(r,c));
            end
            if isnan(chScoreMap(r,c))
                txtcol = 'w';
            else
                txtcol = 'k';
            end
            text(ax2, c, r, txt, 'HorizontalAlignment', 'center', 'Color', txtcol, 'FontSize', 8, 'PickableParts', 'none');
        end
    end
end

%wait for a click, confirm it, repeat until happy (or the figure gets closed)
ishsokay = 0;
while ishsokay == 0
    if ~isvalid(fig) %map got closed while the confirm box was up
        return
    end
    fig.UserData = [];
    uiwait(fig); %chanMapClick drops the clicked channel in fig.UserData and resumes
    if ~isvalid(fig) || isempty(fig.UserData) %closed without picking
        return
    end
    pick = fig.UserData;
    answer = nonModalConfirm(sprintf('Use channel %d?', pick), 'Confirm channel');
    if strcmp(answer, 'Yes')
        chan = pick;
        ishsokay = 1;
    end
end
close(fig);
end

function chanMapClick(img, pmap, chScoreMap, fig)
% ButtonDownFcn for the channel map: turns the click into a (row, shank)
% cell and hands back that cell's raw channel ID via fig.UserData.
% Removed/unscored (NaN) cells don't count - the axes title says so.
ax = ancestor(img, 'axes');
pt = ax.CurrentPoint; %x = shank (column), y = ProbeMaps row
c = round(pt(1,1));
r = round(pt(1,2));
if r < 1 || r > size(pmap,1) || c < 1 || c > size(pmap,2)
    return
end
if isnan(chScoreMap(r,c))
    title(ax, sprintf('CH%d is removed/unscored - pick another', pmap(r,c)));
    return
end
fig.UserData = pmap(r,c);
uiresume(fig);
end

function updateTrConditionalMaster(save_directory, master_filename, tr_conditional, animal, full_filename)
% Keeps <save_directory>\<master_filename> in sync (normally
% tr_conditional_master.mat - see master_filename in this file's header
% doc): a single, small index (table, columns Name/Animal/Filename - no
% TrialIdx, that stays in the per-file tables) of which tr_conditional
% names exist and which animal/filename pairs each has been applied to,
% so it can be scanned on its own without opening every per-recording
% result file. The variable stored inside the file is always named
% tr_conditional_master regardless of master_filename.
masterFile = fullfile(save_directory, master_filename);
if isfile(masterFile)
    master_dat = load(masterFile, 'tr_conditional_master');
    tr_conditional_master = master_dat.tr_conditional_master;
else
    tr_conditional_master = table('Size',[0,3], 'VariableTypes',{'string','string','string'}, 'VariableNames',{'Name','Animal','Filename'});
end

alreadyTracked = strcmp(tr_conditional_master.Name, tr_conditional) & ...
    strcmp(tr_conditional_master.Animal, animal) & ...
    strcmp(tr_conditional_master.Filename, full_filename);
if any(alreadyTracked)
    return
end

newRow = table(string(tr_conditional), string(animal), string(full_filename), ...
    'VariableNames', {'Name','Animal','Filename'});
tr_conditional_master = [tr_conditional_master; newRow];

%% save the master index, and log what was added
fprintf('--- tr_conditional_master: added %s / %s / %s ---\n', tr_conditional, animal, full_filename);
save(masterFile, 'tr_conditional_master');

end
